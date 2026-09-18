"""One-time preprocessing of the Laikipia, Kenya soils (Willard et al.).

Willard, L. L., Muñoz-Carpena, R., Venort, T., Gitonga, J., Maltais-Landry,
G., Caylor, K. K. and Palm, C. A. A high-resolution comprehensive database to
evaluate agricultural impacts at edge of savanna in Kenya. Under review,
Scientific Data.

Reads data/willard_raw/SoilsData.csv (the revised long-format release of
2026-09-18, one row per sample x measurement) and writes
data/willard_reference.csv. Both are git-ignored: the data are unpublished.

What is used, per the data's ReadMe:
  * 86 point samples: grid locations 1 and 4 in each of 43 fields
    (smallholder, large farm or adjacent bush) in Laikipia, Nyeri and Meru,
    January-February 2022; undisturbed cores at 15 cm.
  * Retention: 13 points, pF 0 to 4.2 (SoilWaterCharacteristicCurve_final).
    The first apparatus (sandbox, pF 0-2) did not saturate the cores, and
    the authors raised pF 0-2 by one constant per sample so that pF 0 sits
    at 95 % of porosity, leaving the pressure plates (pF 2.15-4.2) as they
    were. That opens a step of median 0.21 between pF 2 and pF 2.15 -- half
    of each curve's whole water loss, in 85 of 86 samples -- and the tool
    then reads the clays as silty clay loam or loamy sand (6 % exact). Only
    pF 0 (saturation, anchored on porosity) and the pressure plates are
    used here (project decision, 2026-09-18); pF 1-2 are dropped. The
    samples were discarded, so the sandbox cannot be rerun. h_kPa =
    10**pF * 0.0980665, pF 0 taken as h = 0; theta is volume %. Blank
    points (never reported) are dropped.
  * Texture: sieve plus hydrometer, USDA limits (sand > 50 um). The class is
    re-derived from the fractions, as for every other source.
  * Ks: none enters the reference (ksat_cmh is empty; project decision,
    2026-09-18). All three conductivities were measured in the field in the
    dry season, when these vertic clays were at or below wilting point and
    cracked: the permeameter's median is 178 cm/h on the drier clay-rich
    samples against 10-14 cm/h elsewhere, and every method rises with clay.
    The reference's Ks is mostly saturated, swollen lab cores. The field
    values are kept, m/s -> cm/h, as ksat_mpd_cmh (modified Philip-Dunne
    permeameter; the value the authors reject, 3LN06324, is left out),
    ksat_saturo_cmh (dual-head infiltrometer) and k_minidisk_cmh
    (unsaturated, 2 cm suction).
  * Bulk density: mean of the KALRO and UF cores of the sample, the value
    the authors use for porosity; a core they reject is left out.
  * Organic carbon: total carbon of each field's composite (location 0),
    carried to both samples of the field; the soils are not calcareous.
    Organic matter = 1.724 x organic carbon.
  * Coordinates: rounded by the authors to 0.01 degree; none for large farms.
"""

import os

import numpy as np
import pandas as pd

import free_m
import swcc_texture as st

RAW = os.path.join("data", "willard_raw", "SoilsData.csv")
OUT = os.path.join("data", "willard_reference.csv")
CM_TO_KPA = 0.0980665
M_S_TO_CM_H = 100.0 * 3600.0
DEPTH_CM = 15.0
OM_PER_OC = 1.724
MIN_POINTS = 5
MAX_RMSE = 0.03
REJECTED = 20            # QAQC_stat_flag: do not use
PF_SANDBOX = (1.0, 2.0)  # the sandbox points left out (pF 0 is kept)


def read():
    d = pd.read_csv(RAW, low_memory=False)
    d["Text"] = d.Value.astype(str).str.strip()      # class labels
    d["Value"] = pd.to_numeric(d.Value, errors="coerce")
    return d[d.QAQC_stat_flag != REJECTED]


def curve(d):
    """Retention points per sample: {SampleID: (h kPa, theta fraction)}."""
    r = d[d.Measurement.str.match(r"^pF[\d.]+_soilmoisture$")].dropna(
        subset=["Value"])
    pf = r.Measurement.str.extract(r"^pF([\d.]+)_")[0].astype(float)
    keep = ~pf.between(*PF_SANDBOX)
    r, pf = r[keep], pf[keep]
    r = r.assign(h=np.where(pf == 0, 0.0, 10.0 ** pf * CM_TO_KPA),
                 theta=r.Value / 100.0).sort_values(["SampleID", "h"])
    return {s: (g.h.to_numpy(), g.theta.to_numpy())
            for s, g in r.groupby("SampleID")}


def measured_points():
    """{layer_id: (h kPa, theta)} for every sample, for tests that refit."""
    return {f"WL_{s}": p for s, p in curve(read()).items()}


def main():
    if not os.path.exists(RAW):
        raise SystemExit(f"{RAW} not found")
    d = read()
    pts = curve(d)
    wide = d[d.SampleType == "point_sample"].pivot_table(
        index="SampleID", columns="Measurement", values="Value", aggfunc="first")
    site = d.drop_duplicates("SampleID").set_index("SampleID")
    oc = (d[(d.Measurement == "TotalCarbon")]
          .drop_duplicates("SiteID").set_index("SiteID").Value)

    rows = []
    skip = dict(points=0, texture=0, fit=0, rmse=0)
    for s, (h, theta) in pts.items():
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
        bd = np.nanmean([w.get("BulkDensity_KALRO", np.nan),
                         w.get("BulkDensity_UF", np.nan)])
        rows.append(dict(
            layer_id=f"WL_{s}",
            # The two samples of a field share a site and a depth, so they
            # are grouped to keep them out of each other's hold-out folds.
            profile_id=f"WL_{info.SiteID}",
            texture_class=st.usda_class(sand, silt, clay),
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=sand, silt=silt, clay=clay,
            ksat_cmh=np.nan,
            ksat_mpd_cmh=w.get("Ksat_MPD", np.nan) * M_S_TO_CM_H,
            ksat_saturo_cmh=w.get("Ksat_Saturo", np.nan) * M_S_TO_CM_H,
            k_minidisk_cmh=w.get("K_MiniDisk", np.nan) * M_S_TO_CM_H,
            depth_cm=DEPTH_CM, lat=info.Latitude, lon=info.Longitude,
            source_db="Willard_Laikipia", oc=c, om=OM_PER_OC * c,
            porosity=w.get("Porosity", np.nan), bd=bd,
            rmse=rmse, n_points=len(h), **fm,
            sample_type="undisturbed",
            sample_type_source="Willard et al.: retention on undisturbed "
                               "cores (KALRO); Ks in situ (MPD)"))

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"samples with a curve: {len(pts)}")
    print(f"reference layers: {len(out)}   (skipped: {skip})")
    print(f"  fields:         {out.profile_id.nunique()}  "
          f"(with coordinates: {out.lat.notna().sum()} samples)")
    print(f"  field Ks kept aside: permeameter {out.ksat_mpd_cmh.notna().sum()} "
          f"(median {out.ksat_mpd_cmh.median():.1f} cm/h), Saturo "
          f"{out.ksat_saturo_cmh.notna().sum()}, mini-disk "
          f"{out.k_minidisk_cmh.notna().sum()}")
    print(f"  with bd:        {out.bd.notna().sum()} (median {out.bd.median():.2f})")
    print(f"  with oc / om:   {out.oc.notna().sum()}  (median om "
          f"{out.om.median():.1f} %)")
    print(f"  median fit RMSE: {out.rmse.median():.4f}")
    lab = (d[d.Measurement == "TextureClass"].drop_duplicates("SampleID")
           .set_index("SampleID").Text.str.lower())
    agree = out.texture_class.to_numpy() == lab.reindex(
        [i[3:] for i in out.layer_id]).to_numpy()
    print(f"  class from fractions matches the published label: "
          f"{agree.sum()}/{len(out)}")
    print(out.texture_class.value_counts().to_string())


if __name__ == "__main__":
    main()
