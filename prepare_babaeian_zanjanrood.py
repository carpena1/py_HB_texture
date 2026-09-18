"""One-time preprocessing of the Zanjanrood watershed soils (Babaeian et al.).

Babaeian, E., Homaee, M., Vereecken, H., Montzka, C., Norouzi, A. A. and van
Genuchten, M. Th. (2015). A comparative study of multiple approaches for
predicting the soil-water retention curve: hyperspectral information vs.
basic soil properties. Soil Science Society of America Journal 79:1043-1058.
doi:10.2136/sssaj2014.09.0355.

Reads data/babaeian_Zanjanrood_raw/Dataset_SoilPhysicalHydraulicProperties_v2.xlsx
(E. Babaeian's own 2010 measurements, shared with this project for
redistribution) and writes data/babaeian_zanjanrood_reference.csv.

  * 173 samples from 101 locations in the Zanjanrood watershed, north-west
    Iran (Calcixerepts, Haploxerepts, Xerorthents), at 0-15 and 15-30 cm
    (depth_cm is the layer mid-depth, 7.5 or 22.5).
  * Retention: 0-100 cm suction on undisturbed cores (6.8 x 5 cm) on a
    hanging water column; 330-15000 cm on disturbed samples in sand box and
    pressure plates. An undisturbed wet end with a disturbed dry end is the
    standard practice this project counts as undisturbed.
  * Ks: constant head on repacked samples (the data's Readme), cm/d ->
    cm/h. It is used as undisturbed, like the curve: project decision
    (2026-09-14), since the repacking appears to have preserved the range --
    predicted from the rest of the reference, these values come out at the
    right level (median 5.9 against 9.7 cm/h, 86 % within a factor of 10).
  * Texture by sedimentation (USDA limits), organic carbon by Walkley-Black,
    bulk density on paraffin-coated clods.
  * Curves are refitted with the tool's own routine; the file's own fits
    (theta_r fixed at 0.01) are not used. Samples with no point wetter than
    330 cm are skipped, since nothing then constrains saturation.
"""

import os

import numpy as np
import pandas as pd

import free_m
import swcc_texture as st

RAW = os.path.join("data", "babaeian_Zanjanrood_raw",
                   "Dataset_SoilPhysicalHydraulicProperties_v2.xlsx")
OUT = os.path.join("data", "babaeian_zanjanrood_reference.csv")
CM_TO_KPA = 0.0980665
MID_DEPTH = {15: 7.5, 30: 22.5}
OM_PER_OC = 1.724
MIN_POINTS = 5
MAX_RMSE = 0.03
WET_KPA = 6.0            # a curve needs one point wetter than this


def read():
    d = pd.read_excel(RAW, header=1)
    return d[d["Location #"].notna()].reset_index(drop=True)


def curve_columns(d):
    """The retention columns (0, 5 cm, ..., 15000 cm) and their h in kPa."""
    cols = list(d.columns[10:21])
    h_cm = np.array([0.0 if c == 0 else float(str(c).split()[0]) for c in cols])
    return cols, h_cm * CM_TO_KPA


def points(row, cols, h):
    """(h kPa, theta fraction) of one sample, missing values dropped."""
    theta = pd.to_numeric(row[cols], errors="coerce").to_numpy(float)
    ok = np.isfinite(theta)
    return h[ok], theta[ok]


def layer_id(row):
    return f"BB_{row['Location #']}"


def measured_points():
    """{layer_id: (h kPa, theta)} for every sample, for tests that refit."""
    d = read()
    cols, h = curve_columns(d)
    return {layer_id(r): points(r, cols, h) for _, r in d.iterrows()}


def main():
    if not os.path.exists(RAW):
        raise SystemExit(f"{RAW} not found")
    d = read()
    cols, hh = curve_columns(d)

    rows = []
    skip = dict(points=0, wet_end=0, texture=0, fit=0, rmse=0)
    for _, r in d.iterrows():
        h, theta = points(r, cols, hh)
        if len(h) < MIN_POINTS:
            skip["points"] += 1
            continue
        if h.min() > WET_KPA:
            skip["wet_end"] += 1
            continue
        sand, silt, clay = (float(r["% Sand"]), float(r["% Silt"]),
                            float(r["% Clay"]))
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

        oc = float(r["OC (%)"])
        rows.append(dict(
            layer_id=layer_id(r),
            profile_id=f"BB_{str(r['Location #']).split('-')[0]}",
            texture_class=st.usda_class(sand, silt, clay),
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=sand, silt=silt, clay=clay,
            ksat_cmh=float(r["Ks (cm/d)"]) / 24.0,
            depth_cm=MID_DEPTH.get(int(r["Depth (cm)"]), np.nan),
            lat=round(float(r.Lat), 4), lon=round(float(r.Long), 4),
            source_db="Babaeian_Zanjan", oc=oc, om=OM_PER_OC * oc,
            porosity=np.nan, bd=float(r["BD (g/cm3)"]),
            rmse=rmse, n_points=len(h), **fm,
            sample_type="undisturbed",
            sample_type_source="Babaeian et al. 2015: 0-100 cm on undisturbed "
                               "cores, 330-15000 cm on disturbed samples; Ks "
                               "on repacked samples, used as undisturbed "
                               "(project decision)"))

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"samples read:     {len(d)}")
    print(f"reference layers: {len(out)}   (skipped: {skip})")
    print(f"  locations:      {out.profile_id.nunique()}")
    print(f"  with ksat:      {out.ksat_cmh.notna().sum()} "
          f"(median {out.ksat_cmh.median():.2f} cm/h, disturbed samples)")
    print(f"  median fit RMSE: {out.rmse.median():.4f}")
    pub = d.set_index(d.apply(layer_id, axis=1)).Texture.str.strip().str.lower()
    agree = out.texture_class.to_numpy() == pub.reindex(out.layer_id).to_numpy()
    print(f"  class from fractions matches the published label: "
          f"{agree.sum()}/{len(out)}")
    # The file's own fits give alpha in 1/cm.
    own = d.set_index(d.apply(layer_id, axis=1)).reindex(out.layer_id)
    ratio = out.alpha_kpa.to_numpy() / (own.Alpha.to_numpy() / CM_TO_KPA)
    print(f"  refit alpha / published alpha: median {np.nanmedian(ratio):.2f}; "
          f"n median {out.n.median():.3f} vs published {own.n.median():.3f}")
    print(out.texture_class.value_counts().to_string())


if __name__ == "__main__":
    main()
