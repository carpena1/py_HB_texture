"""One-time preprocessing of the Laikipia, Kenya soils (Willard et al.).

Data: Willard, L., Muñoz-Carpena, R., Venort, T., Gitonga, J., Maltais-Landry,
G., Caylor, K. and Palm, C. (2026). A multi-domain field dataset to evaluate
agricultural impacts at the edge of a Kenyan savanna. HydroShare,
doi:10.4211/hs.99080b48206e4b2c82e716516d463a61 (CC BY 4.0). Data descriptor:
Willard et al., Scientific Data (in review).

Reads the published soils table, data/willard_raw/hydroshare/
Soils_Field_Observations/SoilsData.csv (long format, one row per sample x
measurement; download the resource's Soils_Field_Observations/ folder there),
and writes data/willard_reference.csv. The data are the project team's own and
the derived table is distributed; the downloaded files stay local.

What is used, per the published ReadMe (Soils_Field_Observations/ReadMe.txt):
  * 86 point samples: grid locations 1 and 4 in each of 43 fields
    (smallholder, large farm or adjacent bush) in Laikipia, Nyeri and Meru,
    January-February 2022, cores at 15 cm.
  * Retention: the published volumetric series,
    pF<pF>_soilmoisture_volumetric_model (% v/v), 13 points from pF 0 to 4.2.
    It converts the measured gravimetric curve to volumetric water with a bulk
    density that follows the water content: normal shrinkage,
    BD(w) = min(1.6, 1 / (1/2.65 + w)), from the KALRO density of the core
    (ReadMe section 7). These are vertic clays; the conversion and its effect
    are discussed in Vertisol_shrinkage/. Points the authors set to missing
    (rising dry-end values, flag 3) are dropped. pF 0 is taken as h = 0.
    Two caveats from the ReadMe stay with the data: the wet end (pF 0-2) of
    18 samples is short of saturation (flag 23), kept as measured; and
    pF 2.85-4.2 were measured on repacked soil (cf. KSSL's dry end).
  * Texture: sieve plus Bouyoucos hydrometer, USDA limits. The 16 samples of
    the 8 fields flagged for possibly incomplete dispersion (flag 25: andic
    bulk density or pH above 7; their clay may be underestimated) are left
    out, so 70 samples remain. The class is re-derived from the fractions.
  * Ks: the SATURO dual-head infiltrometer (Ksat_Saturo, m/s -> cm/h), 26 of
    the 70 samples, is the sample's Ks (project decision, 2026-10-16). The
    modified Philip-Dunne permeameter and the mini-disk infiltrometer are kept
    in the table (ksat_mpd_cmh, k_minidisk_cmh) but stay out of the reference:
    measured in the field in the dry season, when these clays were cracked,
    the permeameter reads crack flow. Values flagged "do not use" (5) are
    dropped.
  * Bulk density: bd_core is the KALRO density of the core, the one the
    published conversion starts from. The bd column, which the tool matches on
    when a user supplies a bulk density, is the modeled density at 33 kPa
    (pF 2.5: volumetric / gravimetric water), as for the other vertisols and
    as KSSL reports them.
  * Organic carbon: total carbon of each field's composite, carried to both
    samples of the field (the soils are not calcareous apart from the
    high-pH fields left out above). Organic matter = 1.724 x organic carbon.
  * Coordinates: published at 0.01 degree for smallholder fields; none for
    large farms.

History: until 2026-10-31 the table was built from the authors' earlier
releases (the 2026-09-18 file with a rescaled sandbox segment, then L.
Willard's revised curves of 2026-10-16). The published series is the final
version of the latter, with the laboratory's wet-mass basis corrected.
"""

import os

import numpy as np
import pandas as pd

import free_m
import swcc_texture as st

RAW = os.path.join("data", "willard_raw", "hydroshare", "Soils_Field_Observations",
                   "SoilsData.csv")
OUT = os.path.join("data", "willard_reference.csv")
CM_TO_KPA = 0.0980665
M_S_TO_CM_H = 100.0 * 3600.0
OM_PER_OC = 1.724
MIN_POINTS = 5
MAX_RMSE = 0.03
DO_NOT_USE = 5           # QAQC_stat_flag: physically impossible value
DISPERSION = 25          # QAQC_flag: texture may be affected by incomplete dispersion
VOLUMETRIC = r"^pF([\d.]+)_soilmoisture_volumetric_model$"


def read():
    """Point-sample rows of the published table, numeric Value, unusable dropped."""
    d = pd.read_csv(RAW, low_memory=False)
    d["Text"] = d.Value.astype(str).str.strip()        # class labels
    d["Value"] = pd.to_numeric(d.Value, errors="coerce")
    return d[d.QAQC_stat_flag != DO_NOT_USE]


def dispersed(d):
    """Sample ids whose texture carries the dispersion flag."""
    t = d[d.Measurement.isin(["Sand", "Silt", "Clay"]) & (d.QAQC_flag == DISPERSION)]
    return set(t.SampleID)


def curve(d):
    """Published volumetric curves: {SampleID: (h kPa, theta, BD at 33 kPa)},
    the dispersed samples left out."""
    p = d[d.SampleType == "point_sample"]
    skip = dispersed(p)
    vol = p[p.Measurement.str.match(VOLUMETRIC)].dropna(subset=["Value"])
    grav = p[p.Measurement == "pF2.5_soilmoisture_gravimetric"].set_index("SampleID").Value
    out = {}
    for s, g in vol.groupby("SampleID"):
        if s in skip:
            continue
        pf = g.Measurement.str.extract(VOLUMETRIC)[0].astype(float).to_numpy()
        o = np.argsort(pf)
        pf, th = pf[o], g.Value.to_numpy()[o] / 100.0
        h = np.where(pf == 0, 0.0, 10.0 ** pf * CM_TO_KPA)
        th33 = g.Value[g.Measurement == "pF2.5_soilmoisture_volumetric_model"]
        bd33 = float(th33.iloc[0] / grav.get(s, np.nan)) if len(th33) else np.nan
        out[s] = (h, th, bd33)
    return out


def measured_points():
    """{layer_id: (h kPa, theta)} for every kept sample, for tests that refit."""
    return {f"WL_{s}": (h, th) for s, (h, th, _) in curve(read()).items()}


def main():
    if not os.path.exists(RAW):
        raise SystemExit(f"{RAW} not found -- download Soils_Field_Observations/ "
                         "from doi:10.4211/hs.99080b48206e4b2c82e716516d463a61")
    d = read()
    pts = curve(d)
    p = d[d.SampleType == "point_sample"]
    wide = p.pivot_table(index="SampleID", columns="Measurement", values="Value",
                         aggfunc="first")
    site = p.drop_duplicates("SampleID").set_index("SampleID")
    oc = (d[(d.Measurement == "TotalCarbon") & (d.SampleType == "field_aggregate")]
          .drop_duplicates("SiteID").set_index("SiteID").Value)

    rows = []
    skip = dict(points=0, texture=0, fit=0, rmse=0)
    for s, (h, theta, bd33) in pts.items():
        if len(h) < MIN_POINTS:
            skip["points"] += 1
            continue
        w = wide.loc[s]
        sand, silt, clay = w.Sand, w.Silt, w.Clay
        tot = sand + silt + clay
        if not np.isfinite(tot) or not (95 <= tot <= 105):
            skip["texture"] += 1
            continue
        sand, silt, clay = 100 * sand / tot, 100 * silt / tot, 100 * clay / tot

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

        info = site.loc[s]
        c = oc.get(info.SiteID, np.nan)
        rows.append(dict(
            layer_id=f"WL_{s}",
            # The two samples of a field share a site and a depth, so they
            # are grouped to keep them out of each other's hold-out folds.
            profile_id=f"WL_{info.SiteID}",
            texture_class=st.usda_class(sand, silt, clay),
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=sand, silt=silt, clay=clay,
            ksat_cmh=w.get("Ksat_Saturo", np.nan) * M_S_TO_CM_H,
            ksat_mpd_cmh=w.get("Ksat_MPD", np.nan) * M_S_TO_CM_H,
            ksat_saturo_cmh=w.get("Ksat_Saturo", np.nan) * M_S_TO_CM_H,
            k_minidisk_cmh=w.get("K_MiniDisk", np.nan) * M_S_TO_CM_H,
            depth_cm=info.SampleDepth_cm, lat=info.Latitude, lon=info.Longitude,
            source_db="Willard_Laikipia", oc=c, om=OM_PER_OC * c,
            porosity=w.get("Porosity", np.nan), bd=bd33,
            bd_core=w.get("BulkDensity_KALRO", np.nan),
            rmse=rmse, n_points=len(h), **fm,
            sample_type="undisturbed",
            sample_type_source="Willard et al. (HydroShare): retention on cores, "
                               "pF 2.85-4.2 on repacked soil; Ks in situ (Saturo)"))

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"samples with a curve: {len(pts)}  (dispersed, left out: {len(dispersed(p))})")
    print(f"reference layers: {len(out)}   (skipped: {skip})")
    print(f"  fields:         {out.profile_id.nunique()}  "
          f"(with coordinates: {out.lat.notna().sum()} samples)")
    print(f"  Ks (Saturo): {out.ksat_cmh.notna().sum()} (median {out.ksat_cmh.median():.1f} cm/h)")
    print(f"  field Ks kept aside: permeameter {out.ksat_mpd_cmh.notna().sum()} "
          f"(median {out.ksat_mpd_cmh.median():.1f} cm/h), mini-disk "
          f"{out.k_minidisk_cmh.notna().sum()}")
    print(f"  with bd:        {out.bd.notna().sum()} (median at 33 kPa {out.bd.median():.2f}; "
          f"core {out.bd_core.median():.2f})")
    print(f"  with oc / om:   {out.oc.notna().sum()}  (median om "
          f"{out.om.median():.1f} %)")
    print(f"  median fit RMSE: {out.rmse.median():.4f}")
    lab = (p[p.Measurement == "TextureClass"].drop_duplicates("SampleID")
           .set_index("SampleID").Text.str.lower())
    agree = out.texture_class.to_numpy() == lab.reindex(
        [i[3:] for i in out.layer_id]).to_numpy()
    print(f"  class from fractions matches the published label: "
          f"{agree.sum()}/{len(out)}")
    print(out.texture_class.value_counts().to_string())


if __name__ == "__main__":
    main()
