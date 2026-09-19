"""One-time preprocessing of the Yellow River Basin soils (Tong et al. 2024).

Tong, Y., Wang, Y., Zhou, J., Guo, X., Wang, T., Xu, Y., Sun, H., Zhang, P.,
Li, Z. and Lauerwald, R. (2024). Dataset of soil hydraulic parameters in the
Yellow River Basin [dataset]. PANGAEA, doi:10.1594/PANGAEA.965004 (CC BY 4.0);
described in Scientific Data (2024).

Reads data/tong_raw/Tong-etal_2024.tab and writes data/tong_reference.csv.

  * 1,035 retention curves from 475 sites, 0.05-5 m deep, 2008-2019, measured
    by centrifuge at ten suctions from 0.01 to 10 bar (1-1000 kPa), with the
    saturated water content measured. The paper took 2,800 undisturbed cores
    (bulk density, Ks) and 2,925 disturbed samples (texture); the retention
    samples are taken to be the undisturbed cores.
  * Only the authors' van Genuchten fits are published (m = 1-1/n, alpha in
    m^-1), not the measured points, and in 170 curves the fitted thetas is
    more than 0.05 away from the measured one. Each curve is therefore
    rebuilt at ten log-spaced suctions from 1 to 1000 kPa from the authors'
    fit, the measured thetas is added at h = 0, and the whole is refitted with
    the tool's own routine; curves where the two disagree fail the usual
    RMSE limit and are dropped.
  * Texture by laser diffraction (Mastersizer 3000), which reads less clay
    than the pipette and hydrometer methods behind the rest of the reference.
    The column headers give the class limits as 0.02 and 0.5 mm, but the
    fractions reproduce the authors' own USDA classes (98 %), so the limits are
    the USDA 0.002 and 0.05 mm. Samples flagged "Error" for texture are dropped.
  * Ks by constant head on undisturbed cores, cm/min -> cm/h.
"""

import os

import numpy as np
import pandas as pd

import free_m
import swcc_texture as st

RAW = os.path.join("data", "tong_raw", "Tong-etal_2024.tab")
OUT = os.path.join("data", "tong_reference.csv")
M_HEAD_TO_KPA = 9.80665
H_KPA = np.logspace(0.0, 3.0, 10)          # 0.01-10 bar
MAX_RMSE = 0.03
COLS = ["event", "site", "sample", "lon", "lat", "elev", "landuse", "year",
        "depth_m", "clay", "silt", "sand", "psa_method", "psa_quality",
        "class1", "class2", "bd", "bd_method", "ks_cmmin", "ks_method",
        "wrc_method", "ths_meas", "ths", "thr", "alpha_m", "n", "m", "fc",
        "pwp", "r2", "ths_re", "re_range"]


def layer_id(sample):
    return f"TG_{sample}"


def read():
    with open(RAW, encoding="utf-8") as f:
        skip = next(i for i, line in enumerate(f) if line.startswith("*/")) + 1
    d = pd.read_csv(RAW, sep="\t", skiprows=skip)
    d.columns = COLS
    return d[d.n.notna()].reset_index(drop=True)


def measured_points():
    """{layer_id: (h kPa, theta)}: the authors' fit at 1-1000 kPa plus the
    measured saturated water content (the measured points are not published)."""
    out = {}
    for r in read().itertuples():
        th = st.vg_theta(H_KPA, r.thr, r.ths, r.alpha_m / M_HEAD_TO_KPA, r.n)
        out[layer_id(r.sample)] = (np.r_[0.0, H_KPA], np.r_[r.ths_meas, th])
    return out


def main():
    if not os.path.exists(RAW):
        raise SystemExit(f"{RAW} not found")
    d = read()
    pts = measured_points()
    rows = []
    skip = dict(texture=0, fit=0, rmse=0)
    for r in d.itertuples():
        tot = r.sand + r.silt + r.clay
        if r.psa_quality == "Error" or not (95 <= tot <= 105):
            skip["texture"] += 1
            continue
        sand, silt, clay = (100 * v / tot for v in (r.sand, r.silt, r.clay))
        h, theta = pts[layer_id(r.sample)]
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
        rows.append(dict(
            layer_id=layer_id(r.sample), profile_id=f"TG_{r.site}",
            texture_class=st.usda_class(sand, silt, clay),
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=sand, silt=silt, clay=clay,
            ksat_cmh=r.ks_cmmin * 60.0, depth_cm=r.depth_m * 100.0,
            lat=r.lat, lon=r.lon, source_db="Tong_YellowRiver",
            oc=np.nan, om=np.nan, porosity=np.nan, bd=r.bd,
            rmse=rmse, n_points=len(h), **fm,
            sample_type="undisturbed",
            sample_type_source="Tong et al. (2024): undisturbed cores for "
                               "Ks and bulk density; retention taken to be "
                               "on the same cores"))
    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"wrote {OUT}: {len(out)} of {len(d)} curves; dropped {skip}")
    print(out.texture_class.value_counts().to_string())


if __name__ == "__main__":
    main()
