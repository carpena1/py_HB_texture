"""One-time preprocessing of the Fujii et al. soils (Japan and Southeast Asia).

Data: Fujii, K. et al. (2023). Soil water retention curves and soil
physicochemical properties. Dryad, doi:10.5061/dryad.76hdr7svc (CC0); primary
article in Vadose Zone Journal, doi:10.1002/vzj2.70036. The Dryad file
download sits behind a bot check, so the zip is downloaded by hand.

Reads data/fujii_raw/Soil_water_project_tables.xlsx and writes
data/fujii_reference.csv. Under test as an independent set
(verify_external.py fujii): it is not in any reference set the tool uses.

What is used, per the Dryad methods text:
  * 59 horizons of the authors' own sites: Japan 18 (7 sites, several on
    volcanic ash), Indonesia 27 (East Kalimantan), Thailand 12, Canada 2.
    The 10 Czech horizons (Kodesova et al. 2011, compiled from the literature)
    have only 4 points up to 1500 kPa, below MIN_POINTS, so they drop out.
  * Retention on undisturbed 0.1 L cores: saturation and 0.98, 3.1 kPa (sand
    column), 9.8, 31, 98 kPa (pressure plate), pF 4.2 = 1554 kPa
    (centrifuge). The pF 5.5 (air-dry) and pF 7.0 (oven-dry, zero by
    definition) points are dropped: the reference stops at 1500 kPa.
  * Texture by pipette, USDA limits (sand 0.05-2 mm), on the sieved soil.
    The property sheet names some horizons differently from the retention
    sheet, so the two are joined on site and depth.
  * Ks by constant or falling head on the cores. The methods text says m/day
    but the table header says cm/s; the values (median 0.0064, about 23 cm/h
    as cm/s, in forest topsoils with bulk density 0.3-1.0) fit cm/s, so cm/s
    is used.
  * Coordinates from the site table. The East Kalimantan sites MF, AF and
    IG (6 horizons) are not in it, so they have none.
  * Total C by dry combustion (oc, %). Andic per Soil Taxonomy from the
    oxalate extracts: Alo + 1/2 Feo >= 2 % and bulk density <= 0.90.
"""

import os
import sys

# Run from the repository root; the tool's modules sit one folder up.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import re

import numpy as np
import pandas as pd

import free_m
import swcc_texture as st

RAW = os.path.join("data", "fujii_raw", "Soil_water_project_tables.xlsx")
OUT = os.path.join("data", "fujii_reference.csv")
MIN_POINTS = 5
MAX_RMSE = 0.03
# Retention sheet: column index -> suction (kPa); pF 5.5 and 7.0 left out.
H_COL = {7: 0.0, 8: 0.980392, 9: 3.100272, 10: 9.803922, 11: 31.002722,
         12: 98.039216, 13: 1553.816855}
CM_S_TO_CM_H = 3600.0
# Site codes the property sheet spells differently from the retention sheet.
SITE_ALIAS = {"RP": "RPf", "BS": "BSf", "Saskatchewan": "SK06"}


def layer_id(site, hor, depth):
    return f"FJ_{site}_{hor}_{depth}".replace(" ", "")


def _deg(s):
    """'N40˚00΄' -> 40.0; the sign from the hemisphere letter."""
    m = re.match(r"([NSEW])\s*(\d+)\D+(\d+)", s.strip())
    hemi, d, mnt = m.group(1), int(m.group(2)), int(m.group(3))
    if d > 180:              # 'E1135°10'' in the site table
        d = int(str(d)[1:])
    return (-1 if hemi in "SW" else 1) * (d + mnt / 60.0)


def sites():
    """{site code: (lat, lon, parent material)} from the site table."""
    t = pd.read_excel(RAW, "Site information", header=None).iloc[2:]
    t[0] = t[0].ffill()
    out = {}
    for _, r in t[t[2].notna()].iterrows():
        code = str(r[2]).split()[0]
        lat_s, lon_s = str(r[3]).split(",")
        lat, lon = _deg(lat_s), _deg(lon_s)
        if r[0] == "Japan":  # two Japanese sites are printed with W
            lon = abs(lon)
        out[code] = (lat, lon, str(r[8]))
    return out


def read():
    """One row per horizon: retention and property sheets joined."""
    w = pd.read_excel(RAW, "Soil water retention", header=None).iloc[5:]
    w[[0, 2]] = w[[0, 2]].ffill()
    w = w[w[3].notna()]
    w = pd.DataFrame({
        "country": w[0], "site": w[2].astype(str).str.split().str[0],
        "hor": w[3].astype(str).str.strip(),
        "depth": w[4].astype(str).str.strip().str.rstrip("+"),
        **{f"th_{c}": pd.to_numeric(w[c], errors="coerce") for c in H_COL},
        "bd": pd.to_numeric(w[17], errors="coerce"),
        "ks_cms": pd.to_numeric(w[20], errors="coerce")})
    p = pd.read_excel(RAW, "Soil property", header=None).iloc[4:]
    p = p[p[2].notna()].copy()               # country lines carry no horizon
    p[0] = p[0].ffill()
    p = pd.DataFrame({
        "site": p[0].astype(str).str.split().str[0].str.rstrip(",").replace(SITE_ALIAS),
        "hor_prop": p[2].astype(str).str.strip(),
        "depth": p[4].astype(str).str.strip().str.rstrip("+"),
        "tc_gkg": pd.to_numeric(p[8], errors="coerce"),
        "sand": pd.to_numeric(p[24], errors="coerce"),
        "silt": pd.to_numeric(p[26], errors="coerce"),
        "clay": pd.to_numeric(p[28], errors="coerce"),
        "texture_lab": p[30], "feo_gkg": pd.to_numeric(p[38], errors="coerce"),
        "alo_gkg": pd.to_numeric(p[40], errors="coerce"),
        "allophane_pct": pd.to_numeric(p[48], errors="coerce"),
        "data_source": p[55]})
    # The two sheets name some horizons differently (A1, A2, BA against A, AB,
    # B at the same depths), so they are joined on site and depth.
    d = w.merge(p, on=["site", "depth"], how="left", validate="one_to_one")
    d["layer_id"] = [layer_id(s, h, z) for s, h, z in zip(d.site, d.hor, d.depth)]
    return d.reset_index(drop=True)


def points(r):
    h = np.array(list(H_COL.values()))
    th = np.array([r[f"th_{c}"] for c in H_COL], dtype=float)
    ok = np.isfinite(th)
    return h[ok], th[ok]


def measured_points():
    """{layer_id: (h kPa, theta)} for every horizon with enough points."""
    out = {}
    for _, r in read().iterrows():
        h, th = points(r)
        if len(h) >= MIN_POINTS:
            out[r.layer_id] = (h, th)
    return out


def main():
    if not os.path.exists(RAW):
        raise SystemExit(f"{RAW} not found -- download the zip from "
                         "https://doi.org/10.5061/dryad.76hdr7svc")
    d = read()
    geo = sites()
    rows = []
    skip = dict(points=0, texture=0, fit=0, rmse=0)
    for _, r in d.iterrows():
        h, theta = points(r)
        if len(h) < MIN_POINTS:
            skip["points"] += 1
            continue
        tot = r.sand + r.silt + r.clay
        if not np.isfinite(tot) or not (95 <= tot <= 105):
            skip["texture"] += 1
            continue
        sand, silt, clay = (100 * r[c] / tot for c in ("sand", "silt", "clay"))
        try:
            popt, _ = st.fit_vg(h, theta)
        except Exception:
            skip["fit"] += 1
            continue
        tr, ts, la, ln1 = popt
        alpha, n = 10.0 ** la, 1.0 + 10.0 ** ln1
        rmse = float(np.sqrt(np.mean(
            (st.vg_theta(h, tr, ts, alpha, n) - theta) ** 2)))
        if not np.isfinite(rmse) or rmse > MAX_RMSE or n <= 1.0:
            skip["rmse"] += 1
            continue
        fm = free_m.fit_free_m(h, theta,
                               fallback=free_m.mualem_fallback(tr, ts, alpha, n))
        top, bot = (float(x) for x in r.depth.split("-"))
        lat, lon, parent = geo.get(r.site, (np.nan, np.nan, ""))
        alo_feo = (r.alo_gkg + 0.5 * r.feo_gkg) / 10.0
        andic = ("yes" if alo_feo >= 2.0 and r.bd <= 0.90 else
                 "no" if np.isfinite(alo_feo) and np.isfinite(r.bd) else np.nan)
        rows.append(dict(
            layer_id=r.layer_id, profile_id=f"FJ_{r.site}",
            texture_class=st.usda_class(sand, silt, clay),
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=sand, silt=silt, clay=clay,
            ksat_cmh=r.ks_cms * CM_S_TO_CM_H if r.ks_cms > 0 else np.nan,
            depth_cm=0.5 * (top + bot), lat=lat, lon=lon, country=r.country,
            source_db="Fujii_Dryad", oc=r.tc_gkg / 10.0, bd=r.bd,
            alo_feo_pct=alo_feo, allophane_pct=r.allophane_pct,
            # not 'andic': load_reference_df merges that from the volcanic tables
            andic_st=andic,
            parent_material=parent, texture_lab=r.texture_lab,
            rmse=rmse, n_points=len(h), **fm,
            sample_type="undisturbed",
            sample_type_source="Dryad methods: retention on undisturbed 0.1 L cores"))
    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"horizons: {len(d)}; wrote {OUT}: {len(out)} from "
          f"{out.profile_id.nunique()} sites (skipped: {skip})")
    print(f"  by country: {out.country.value_counts().to_dict()}")
    print(f"  Ks: {out.ksat_cmh.notna().sum()} (median {out.ksat_cmh.median():.1f} cm/h); "
          f"andic: {(out.andic_st == 'yes').sum()}; median fit RMSE {out.rmse.median():.4f}")
    lab = out.texture_lab.astype(str).str.lower()
    print(f"  class from fractions matches the file's class: "
          f"{(lab == out.texture_class).sum()}/{len(out)}")
    print(out.texture_class.value_counts().to_string())


if __name__ == "__main__":
    main()
