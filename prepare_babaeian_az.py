"""One-time preprocessing of E. Babaeian's Arizona soils.

Reads data/babaeian_AZ_raw/AZ_SoilData_Basic_SWC_vanGenuchten_parameters.xlsx
(E. Babaeian's own measurements, shared with this project) and writes
data/babaeian_az_reference.csv.

  * 21 samples, AZ1-AZ20 with AZ4 as A and B, from Arizona, USA; sand to clay.
    No coordinates or depths are given (lat, lon, depth_cm left empty).
    AZ4A and AZ4B share one profile_id.
  * Retention: 0.01 cm (saturation) and 7.4-800 cm, then a dry end of
    repeated readings from about 1,000 to 2,900,000 cm. Only points up to
    1500 kPa are used: the reference is measured and described over that
    range, and the dry end (up to 18 of a sample's 27 points) would otherwise
    dominate the fit. AZ18 fits at RMSE 0.008 that way against 0.032 with
    every point.
  * The 500 cm reading of AZ15, AZ18 and AZ19 falls below the 800 cm one
    (0.215, 0.304, 0.219); it is replaced by interpolating, linear in
    log10(h), between the 150 and 800 cm readings (project decision,
    2026-09-17). With that all 21 samples fit.
  * Samples are undisturbed (confirmed to the project, 2026-09-17).
  * Ks in cm/d -> cm/h; organic matter as given, oc = om / 1.724.
  * Curves are refitted with the tool's own routine; the file's own fits are
    not used.
"""

import os

import numpy as np
import pandas as pd

import free_m
import swcc_texture as st

RAW = os.path.join("data", "babaeian_AZ_raw",
                   "AZ_SoilData_Basic_SWC_vanGenuchten_parameters.xlsx")
OUT = os.path.join("data", "babaeian_az_reference.csv")
CM_TO_KPA = 0.0980665
OM_PER_OC = 1.724
MIN_POINTS = 5
MAX_RMSE = 0.03
MAX_KPA = 1500.0         # the range the reference covers
# Samples whose 500 cm reading is out of order, interpolated from 150 and 800 cm.
FIX_500CM = ("AZ15", "AZ18", "AZ19")


def layer_id(sample):
    return f"BA_{sample}"


def fix_500cm(h_cm, theta):
    """Replace the 500 cm reading by log10(h)-linear interpolation between
    the 150 and 800 cm readings."""
    i150, i500, i800 = (int(np.argmin(np.abs(h_cm - v)))
                        for v in (150.0, 500.1, 800.0))
    w = ((np.log10(h_cm[i500]) - np.log10(h_cm[i150]))
         / (np.log10(h_cm[i800]) - np.log10(h_cm[i150])))
    theta = theta.copy()
    theta[i500] = theta[i150] + w * (theta[i800] - theta[i150])
    return theta


def measured_points(max_kpa=MAX_KPA):
    """{layer_id: (h kPa, theta)} for every sample, for tests that refit;
    max_kpa=None keeps the dry end."""
    d = pd.read_excel(RAW, sheet_name=0, header=None)
    out = {}
    for c in range(1, d.shape[1], 2):
        sample = str(d.iloc[0, c]).strip()
        h = pd.to_numeric(d.iloc[2:, c], errors="coerce").to_numpy(float)
        th = pd.to_numeric(d.iloc[2:, c + 1], errors="coerce").to_numpy(float)
        ok = np.isfinite(h) & np.isfinite(th)
        h, th = h[ok], th[ok]
        if sample in FIX_500CM:
            th = fix_500cm(h, th)
        if max_kpa is not None:
            keep = h * CM_TO_KPA <= max_kpa
            h, th = h[keep], th[keep]
        out[layer_id(sample)] = (h * CM_TO_KPA, th)
    return out


def properties():
    d = pd.read_excel(RAW, sheet_name=1, header=1)
    d = d[d["Sample_ID"].astype(str).str.startswith("AZ")]
    return d.reset_index(drop=True)


def main():
    if not os.path.exists(RAW):
        raise SystemExit(f"{RAW} not found")
    pts = measured_points()
    props = properties()

    rows = []
    skip = dict(points=0, texture=0, fit=0, rmse=0)
    for _, r in props.iterrows():
        sample = str(r["Sample_ID"]).strip()
        h, theta = pts[layer_id(sample)]
        if len(h) < MIN_POINTS:
            skip["points"] += 1
            continue
        sand, silt, clay = (float(r["Sand (%)"]), float(r["Silt (%)"]),
                            float(r["Clay (%)"]))
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

        om = float(r["Organic matter (%w)"])
        rows.append(dict(
            layer_id=layer_id(sample),
            profile_id=f"BA_{sample.rstrip('AB')}",
            texture_class=st.usda_class(sand, silt, clay),
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=sand, silt=silt, clay=clay,
            ksat_cmh=float(r["Ks (cm/d)"]) / 24.0,
            depth_cm=np.nan, lat=np.nan, lon=np.nan,
            source_db="Babaeian_Arizona", oc=om / OM_PER_OC, om=om,
            porosity=np.nan, bd=float(r["Bulk density (g/cm3)"]),
            rmse=rmse, n_points=len(h), **fm,
            sample_type="undisturbed",
            sample_type_source="undisturbed samples (confirmed to the "
                               "project, 2026-09-17)"))

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"samples read:     {len(props)}")
    print(f"reference layers: {len(out)}   (skipped: {skip})")
    print(f"  with ksat:      {out.ksat_cmh.notna().sum()} "
          f"(median {out.ksat_cmh.median():.2f} cm/h)")
    print(f"  median fit RMSE: {out.rmse.median():.4f}")
    pub = props.set_index(props.Sample_ID.map(layer_id))[
        "USDA Textural Class"].str.strip().str.lower()
    agree = out.texture_class.to_numpy() == pub.reindex(out.layer_id).to_numpy()
    print(f"  class from fractions matches the published label: "
          f"{agree.sum()}/{len(out)}")
    print(out.texture_class.value_counts().to_string())


if __name__ == "__main__":
    main()
