"""Volcanic parent material and andic properties for every reference layer.

Writes two tables keyed on layer_id, which load_reference_df() joins when
they exist:

  data/volcanic_flags.csv             layers of the distributed tables
  data/volcanic_flags_restricted.csv  layers of restricted or unpublished
                                      tables (EU-HYDI, Laikipia, Canary,
                                      Arizona); git-ignored

Two attributes, kept apart because they differ: the Laikipia soils sit on
Mount Kenya's volcanics but show no andic properties.

volcanic -- volcanic parent material
  known     the source says so: KSSL Andisols, andic or vitrandic subgroups
            and layers meeting the USDA andic lab criteria; Africa Soil
            Profiles (AfSP) profiles on volcanic rock, ash or tuff, or
            classified Andosols; the Canary and Laikipia sets; Romano's
            Campania sites (EU-HYDI report EUR 26053, chapter 14: "andic
            features ... affected by ... Vesuvio or Vulture"); the Russian
            register's Kamchatka profiles (a peninsula mantled by volcanic
            ash; all within 75 km of an active volcano). The register's
            per-profile classes would be better, but its only download
            (egrpr.esoil.ru soil_data.xls) is truncated at the source.
  probable  SoilGrids 2.0 Andosols probability >= 20 % or Andosols most
            probable, or a Holocene volcano within 25 km.
  possible  any Holocene or Pleistocene volcano within 100 km and SoilGrids
            Andosols probability >= 5 %.
  unlikely  none of the above, or the source rules it out: a KSSL taxonomy
            that is not andic; an AfSP profile with a non-volcanic parent
            material or class; the EU-HYDI Sicily chapter (15: no andic
            soils); Romano's Upper Alento sites (flysch clays).
  unknown   no coordinates and no source information.

andic -- andic soil properties
  yes       KSSL andic taxonomy or lab criteria; AfSP Andosols; Canary.
  likely    Romano's Campania sites (the source's own description);
            Kamchatka profiles whose median bulk density is <= 0.90, the
            andic limit.
  no        a KSSL or AfSP classification that is not andic; Laikipia
            (bulk density ~1.05, water at 1500 kPa 0.35 x clay).
  unknown   everything else.

Checked against the KSSL labels (2026-09-17): SoilGrids is specific but
insensitive (p >= 20 %: 15 % of andic layers found, 99 % of others
excluded), and distance alone flags arid volcanic fields that never formed
andic soils. Hence the tiers.

Inputs (all git-ignored, see .gitignore): data/kssl_raw (NCSS Lab Data
Mart), data/volcanoes_raw (Smithsonian GVP Holocene and Pleistocene lists,
VOTW 5.4.0), data/soil_registries_raw/AF-AfSP1.2.zip (ISRIC), and the
SoilGrids 2.0 WRB layers read over the network once and cached in
data/volcanic_raw/.

Usage:  python prepare_volcanic.py
"""

import glob
import io
import os
import re
import sqlite3
import struct
import sys
import time
import xml.etree.ElementTree as ET
import zipfile

import numpy as np
import pandas as pd

import swcc_texture as st

CACHE = os.path.join("data", "volcanic_raw")
OUT_PUBLIC = os.path.join("data", "volcanic_flags.csv")
OUT_RESTRICTED = os.path.join("data", "volcanic_flags_restricted.csv")
KSSL_DB = os.path.join("data", "kssl_raw", "NCSSLabDataMartSQLite.sqlite3")
AFSP_ZIP = os.path.join("data", "soil_registries_raw", "AF-AfSP1.2.zip")
SOILGRIDS = "/vsicurl/https://files.isric.org/soilgrids/latest/data/wrb/"
WRB = ("Acrisols Albeluvisols Alisols Andosols Arenosols Calcisols Cambisols "
       "Chernozems Cryosols Durisols Ferralsols Fluvisols Gleysols Gypsisols "
       "Histosols Kastanozems Leptosols Lixisols Luvisols Nitisols Phaeozems "
       "Planosols Plinthosols Podzols Regosols Solonchaks Solonetz Stagnosols "
       "Umbrisols Vertisols").split()
EARTH_KM = 6371.0

# Source-document rules (lat, lon, radius km).
ROMANO_CAMPANIA = [(40.94, 14.37, 10), (40.93, 14.20, 10), (40.98, 14.26, 10),
                   (40.99, 14.19, 10), (40.92, 14.71, 10)]   # Acerra, Giugliano,
                                                            # Succivo, Lusciano,
                                                            # Monteforte Irpino
ROMANO_ALENTO = [(40.38, 15.18, 25)]
KAMCHATKA = dict(lat=(50.8, 60.5), lon=(155.5, 163.5))
PUBLIC_SOURCES = {"gshp", "kssl", "hohenbrink", "babaeian_zanjanrood",
                  "unsoda", "sdb"}


def layers():
    frames = []
    for name in ("gshp", "kssl", "hohenbrink", "babaeian_zanjanrood", "unsoda",
                 "sdb", "euhydi", "willard", "armas", "babaeian_az"):
        path = os.path.join(st.DATA_DIR, f"{name}_reference.csv")
        if os.path.exists(path):
            frames.append(pd.read_csv(path, low_memory=False).assign(table=name))
    d = pd.concat(frames, ignore_index=True)
    d["source"] = d.source_db.fillna(d.layer_id.str.split("_").str[0])
    d.loc[d.lat.abs() > 90, "lat"] = np.nan
    return d.drop_duplicates("layer_id").reset_index(drop=True)


def kssl_taxonomy():
    """Per KSSL layer: classified?, andic by taxonomy or by the lab criteria."""
    q = """SELECT l.layer_key, c.samp_taxorder, c.corr_taxorder, c.SSL_taxorder,
                  c.samp_taxsubgrp, c.corr_taxsubgrp, c.SSL_taxsubgrp,
                  ch.aluminum_plus_half_iron_oxalate AS alfe,
                  ch.new_zealand_phosphorus_retent AS pret,
                  p.bulk_density_third_bar AS bd33
           FROM lab_layer l
           LEFT JOIN lab_combine_nasis_ncss c ON c.pedon_key = l.pedon_key
           LEFT JOIN lab_chemical_properties_vw ch ON ch.layer_key = l.layer_key
           LEFT JOIN lab_physical_properties_vw p ON p.layer_key = l.layer_key"""
    with sqlite3.connect(KSSL_DB) as con:
        k = pd.read_sql(q, con)
    k["layer_id"] = "KSSL_" + k.layer_key.astype(str)
    orders = k[["samp_taxorder", "corr_taxorder", "SSL_taxorder"]].apply(
        lambda s: s.str.lower())
    subs = k[["samp_taxsubgrp", "corr_taxsubgrp", "SSL_taxsubgrp"]].apply(
        lambda s: s.str.lower())
    k["classified"] = orders.notna().any(axis=1)
    k["andisol"] = (orders == "andisols").any(axis=1)
    k["andic_sub"] = subs.apply(lambda s: s.str.contains(r"andic|vitrand",
                                                         na=False)).any(axis=1)
    k["lab_andic"] = ((pd.to_numeric(k.alfe, errors="coerce") >= 2.0)
                      & (pd.to_numeric(k.pret, errors="coerce") >= 85)
                      & (pd.to_numeric(k.bd33, errors="coerce") <= 0.90))
    k["andic_tax"] = k.andisol | k.andic_sub | k.lab_andic
    return k.drop_duplicates("layer_id")[["layer_id", "classified", "andic_tax",
                                          "andisol", "andic_sub", "lab_andic"]]


def soilgrids(pts):
    """Andosols probability and most probable WRB group at each point; cached."""
    path = os.path.join(CACHE, "soilgrids_points.csv")
    have = pd.read_csv(path) if os.path.exists(path) else pd.DataFrame(
        columns=["la", "lo", "p_andosols", "wrb_code"])
    need = pts.merge(have[["la", "lo"]], how="left", indicator=True)
    need = need[need._merge == "left_only"][["la", "lo"]].reset_index(drop=True)
    if len(need):
        import rasterio
        from rasterio.warp import transform
        print(f"SoilGrids: reading {len(need)} new points over the network",
              flush=True)
        env = dict(GDAL_HTTP_MULTIRANGE="YES",
                   GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES", VSI_CACHE="TRUE",
                   VSI_CACHE_SIZE=str(512 * 1024 * 1024),
                   GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                   GDAL_HTTP_MAX_RETRY="5", GDAL_HTTP_RETRY_DELAY="2")
        order = np.lexsort((need.lo.values, need.la.values))
        with rasterio.Env(**env):
            for name, col in (("Andosols.vrt", "p_andosols"),
                              ("MostProbable.vrt", "wrb_code")):
                t0 = time.time()
                with rasterio.open(SOILGRIDS + name) as src:
                    xs, ys = transform("EPSG:4326", src.crs,
                                       need.lo.values[order].tolist(),
                                       need.la.values[order].tolist())
                    vals = np.empty(len(need))
                    for i, v in enumerate(src.sample(list(zip(xs, ys)))):
                        vals[order[i]] = v[0]
                        if i % 500 == 0:
                            print(f"  {name} {i}/{len(need)} "
                                  f"({time.time() - t0:.0f} s)", flush=True)
                need[col] = vals
        have = pd.concat([have, need], ignore_index=True)
        os.makedirs(CACHE, exist_ok=True)
        have.to_csv(path, index=False)
    have = have.copy()
    have.loc[have.p_andosols > 100, "p_andosols"] = np.nan        # nodata
    have["wrb"] = [WRB[int(c)] if np.isfinite(c) and 0 <= c < len(WRB)
                   else None for c in have.wrb_code.astype(float)]
    return have


def gvp_volcanoes():
    """Smithsonian GVP lists (SpreadsheetML; stray '<' in text is escaped)."""
    ns = {"ss": "urn:schemas-microsoft-com:office:spreadsheet"}
    out = []
    for p in sorted(glob.glob(os.path.join("data", "volcanoes_raw",
                                           "GVP_Volcano_List_*.xls"))):
        txt = open(p, "rb").read().decode("utf-8", "replace")
        txt = re.sub(r"<(?![/A-Za-z?!])", "&lt;", txt)
        root = ET.fromstring(txt)
        rows = []
        for r in root.find("ss:Worksheet", ns).iter(f"{{{ns['ss']}}}Row"):
            vals = []
            for c in r.findall("ss:Cell", ns):
                idx = c.get(f"{{{ns['ss']}}}Index")
                while idx and len(vals) < int(idx) - 1:
                    vals.append(None)
                d = c.find("ss:Data", ns)
                vals.append(d.text if d is not None else None)
            rows.append(vals)
        hdr = rows[1]
        df = pd.DataFrame([r + [None] * (len(hdr) - len(r)) for r in rows[2:]],
                          columns=hdr)
        df["holocene"] = "Holocene" in p
        out.append(df)
    v = pd.concat(out, ignore_index=True)
    for c in ("Latitude", "Longitude"):
        v[c] = pd.to_numeric(v[c], errors="coerce")
    return v.dropna(subset=["Latitude", "Longitude"]).reset_index(drop=True)


def nearest_km(la, lo, pla, plo):
    from sklearn.neighbors import BallTree
    tree = BallTree(np.radians(np.c_[pla, plo]), metric="haversine")
    d, i = tree.query(np.radians(np.c_[la, lo]), k=1)
    return d[:, 0] * EARTH_KM, i[:, 0]


def afsp_profiles():
    """Africa Soil Profiles 1.2: coordinates, WRB/FAO/USDA class, parent
    material, from the zipped dBase table (read without extra libraries)."""
    path = os.path.join(CACHE, "afsp_profiles.csv")
    if os.path.exists(path):
        return pd.read_csv(path, low_memory=False)
    with zipfile.ZipFile(AFSP_ZIP) as z:
        raw = z.read("AfSP012Qry_ISRIC/GIS_Dbf/AfSP012Qry_Profiles.dbf")
    nrec, hlen, rlen = struct.unpack("<I H H", raw[4:12])
    fields, pos = {}, 1
    for k in range(32, hlen - 1, 32):
        name = raw[k:k + 11].split(b"\0")[0].decode("latin-1")
        fields[name] = (pos, raw[k + 16])
        pos += raw[k + 16]
    want = ["ProfileID", "X_LonDD", "Y_LatDD", "WRB06", "FAO88", "USDA",
            "ParMat", "Litholo"]
    rows = []
    for r in range(nrec):
        rec = raw[hlen + r * rlen: hlen + (r + 1) * rlen]
        rows.append([rec[fields[w][0]:fields[w][0] + fields[w][1]]
                     .decode("latin-1").strip() for w in want])
    p = pd.DataFrame(rows, columns=want).replace({"NA": np.nan, "": np.nan})
    for c in ("X_LonDD", "Y_LatDD"):
        p[c] = pd.to_numeric(p[c], errors="coerce")
    os.makedirs(CACHE, exist_ok=True)
    p.to_csv(path, index=False)
    return p


def within(df, sites):
    m = np.zeros(len(df), bool)
    for la, lo, r in sites:
        m |= (np.hypot((df.lat - la) * 111.0,
                       (df.lon - lo) * 111.0 * np.cos(np.radians(la))) <= r
              ).to_numpy()
    return m


def main():
    d = layers()
    print(f"layers: {len(d)} from {d.table.nunique()} tables")
    geo = d.lat.notna() & d.lon.notna()
    d["la"], d["lo"] = d.lat.round(2), d.lon.round(2)

    sg = soilgrids(d.loc[geo, ["la", "lo"]].drop_duplicates())
    d = d.merge(sg[["la", "lo", "p_andosols", "wrb"]], on=["la", "lo"], how="left")

    v = gvp_volcanoes()
    for tag, sel in (("any", np.ones(len(v), bool)), ("holocene", v.holocene)):
        km, i = nearest_km(d.loc[geo, "lat"], d.loc[geo, "lon"],
                           v.Latitude[sel], v.Longitude[sel])
        d.loc[geo, f"km_{tag}"] = km
        d.loc[geo, f"volcano_{tag}"] = v["Volcano Name"][sel].to_numpy()[i]

    k = kssl_taxonomy()
    d = d.merge(k, on="layer_id", how="left")
    kssl_cls = d.classified.eq(True)
    kssl_andic = d.andic_tax.eq(True)

    af = afsp_profiles().dropna(subset=["X_LonDD", "Y_LatDD"])
    isaf = (d.source == "AfSPDB") & geo
    km, i = nearest_km(d.loc[isaf, "lat"], d.loc[isaf, "lon"],
                       af.Y_LatDD, af.X_LonDD)
    m = af.iloc[i].reset_index(drop=True)
    close = km <= 1.0
    af_andosol = (m.WRB06.str.contains("Andosol", case=False, na=False)
                  | m.FAO88.str.contains("Andosol", case=False, na=False)
                  | m.USDA.str.contains(r"\bAndisol|and(?:s|epts?|ists?)\b",
                                        case=False, regex=True, na=False))
    af_volc = (m.ParMat.str.match(r"^V", na=False)
               | m.Litholo.str.match(r"^V", na=False))
    af_class = m[["WRB06", "FAO88", "USDA", "ParMat", "Litholo"]].notna().any(axis=1)
    idx = d.index[isaf]
    d["af_volcanic"], d["af_andosol"], d["af_classified"] = np.nan, np.nan, np.nan
    d.loc[idx, "af_volcanic"] = np.where(close, af_volc | af_andosol, np.nan)
    d.loc[idx, "af_andosol"] = np.where(close, af_andosol, np.nan)
    d.loc[idx, "af_classified"] = np.where(close, af_class, np.nan)

    campania = (d.source == "EUHYDI_Romano") & within(d, ROMANO_CAMPANIA)
    alento = (d.source == "EUHYDI_Romano") & within(d, ROMANO_ALENTO)
    sicily = d.source == "EUHYDI_Iovino"
    kamchatka = ((d.source == "Russia_EGRPR") & d.lat.between(*KAMCHATKA["lat"])
                 & d.lon.between(*KAMCHATKA["lon"]))
    light = d.groupby("profile_id").bd.transform("median") <= 0.90

    known = (kssl_andic | d.af_volcanic.eq(1) | campania
             | kamchatka
             | d.source.isin(["Armas_Canarias", "Willard_Laikipia"]))
    ruled_out = ((kssl_cls & ~kssl_andic) | d.af_volcanic.eq(0) | alento | sicily)
    probable = ((d.p_andosols >= 20) | (d.wrb == "Andosols")
                | (d.km_holocene <= 25))
    possible = (d.km_any <= 100) & (d.p_andosols >= 5)
    d["volcanic"] = np.select(
        [known, ruled_out, probable, possible, ~geo],
        ["known", "unlikely", "probable", "possible", "unknown"], "unlikely")

    d["andic"] = np.select(
        [kssl_andic | d.af_andosol.eq(1) | (d.source == "Armas_Canarias"),
         campania | (kamchatka & light),
         (kssl_cls & ~kssl_andic) | (d.af_classified.eq(1) & d.af_andosol.eq(0))
         | (d.source == "Willard_Laikipia")],
        ["yes", "likely", "no"], "unknown")

    ev = np.full(len(d), "", object)
    rules = [(d.andisol.eq(True), "KSSL Andisol"),
             (d.andic_sub.eq(True), "KSSL andic subgroup"),
             (d.lab_andic.eq(True), "KSSL andic lab criteria"),
             (kssl_cls & ~kssl_andic, "KSSL taxonomy not andic"),
             (d.af_andosol.eq(1), "AfSP Andosol"),
             (d.af_volcanic.eq(1) & d.af_andosol.ne(1), "AfSP volcanic parent material"),
             (d.af_volcanic.eq(0), "AfSP not volcanic"),
             (campania, "EU-HYDI ch.14, Campania"),
             (alento, "Upper Alento flysch"), (sicily, "EU-HYDI ch.15, Sicily"),
             (kamchatka, "Kamchatka volcanic ash"),
             (kamchatka & light, "bulk density <= 0.90"),
             (d.source == "Armas_Canarias", "Canary andic set"),
             (d.source == "Willard_Laikipia", "Mount Kenya volcanics, not andic"),
             ((d.p_andosols >= 20) | (d.wrb == "Andosols"), "SoilGrids Andosols"),
             (d.km_holocene <= 25, "Holocene volcano <= 25 km"),
             ((d.km_any <= 100) & (d.p_andosols >= 5), "volcano <= 100 km + SoilGrids >= 5 %")]
    for mask, text in rules:
        mask = np.asarray(mask, bool)
        ev[mask] = [f"{e}; {text}" if e else text for e in ev[mask]]
    d["volcanic_evidence"] = ev

    cols = ["layer_id", "volcanic", "andic", "volcanic_evidence", "p_andosols",
            "wrb", "km_any", "volcano_any", "km_holocene", "volcano_holocene"]
    pub = d.table.isin(PUBLIC_SOURCES)
    d.loc[pub, cols].round(1).to_csv(OUT_PUBLIC, index=False)
    d.loc[~pub, cols].round(1).to_csv(OUT_RESTRICTED, index=False)
    print(f"wrote {OUT_PUBLIC} ({pub.sum()} layers) and {OUT_RESTRICTED} "
          f"({(~pub).sum()})")
    print(pd.crosstab(d.volcanic, d.andic, margins=True).to_string())
    t = d[d.volcanic.isin(["known", "probable", "possible"])]
    print(pd.crosstab(t.source, t.volcanic).to_string())


if __name__ == "__main__":
    main()
