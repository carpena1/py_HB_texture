"""One-time preprocessing of the Hohenbrink et al. (2023) dataset.

Hohenbrink, T. L., Jackisch, C., Durner, W., Germer, K., Iden, S. C.,
Kreiselmeier, J., Leuther, F., Metzger, J. C., Naseri, M., Peters, A. (2023).
"Soil water retention and hydraulic conductivity measured in a wide saturation
range." Earth System Science Data.

Reads data/Hohenbrink_raw/{BasicProp,MetaData,RetMeas,CondMeas}.csv and writes
data/hohenbrink_reference.csv.

Why this dataset matters here: it is measured by HYPROP evaporation plus WP4
dewpoint, so single curves run from pF -0.03 (effectively saturation) to
pF > 5. Existing reference tables are thin exactly at the wet end -- only
~12 % of GSHP curves have a point at zero suction, and KSSL begins at 6 kPa.

POINT THINNING (important). HYPROP returns ~199 points between pF -0.03 and
~2.6, while the dry limb carries only ~3 WP4 points. Fitting all of them by
least squares weights the wet-to-mid range about 98:2 against the dry end --
which is the part that carries texture, and the part the free m governs. So
points are binned in pF and one median point is kept per bin, giving roughly
even weight per decade of suction.

Units:
    pF        log10(|h| in cm)  ->  h_kPa = 10**pF * 0.0980665
    k         cm/day            ->  ksat_cmh = k / 24
"""

import os

import numpy as np
import pandas as pd

import free_m
import swcc_texture as st

RAW = os.path.join("data", "Hohenbrink_raw")
OUT = os.path.join("data", "hohenbrink_reference.csv")
CM_TO_KPA = 0.0980665
PF_BIN = 0.25            # width of the log-suction bins used for thinning
MIN_POINTS = 6
MAX_RMSE = 0.03
MAX_THETA = 0.95

# BasicProp TexClass_USDA codes -> the names used everywhere else.
TEX = {"Sa": "sand", "LoSa": "loamy sand", "SaLo": "sandy loam",
       "Lo": "loam", "Si": "silt", "SiLo": "silt loam",
       "SaClLo": "sandy clay loam", "ClLo": "clay loam",
       "SiClLo": "silty clay loam", "SaCl": "sandy clay",
       "SiCl": "silty clay", "Cl": "clay"}


def thin(pf, theta):
    """One median point per PF_BIN-wide bin of log suction."""
    d = pd.DataFrame({"pf": pf, "theta": theta})
    d = d[np.isfinite(d.pf) & np.isfinite(d.theta)]
    key = np.floor(d.pf / PF_BIN)
    g = d.groupby(key).median()
    return g.pf.to_numpy(), g.theta.to_numpy()


def main():
    if not os.path.isdir(RAW):
        raise SystemExit(f"{RAW} not found")
    basic = pd.read_csv(os.path.join(RAW, "BasicProp.csv"))
    meta = pd.read_csv(os.path.join(RAW, "MetaData.csv"))
    ret = pd.read_csv(os.path.join(RAW, "RetMeas.csv"))
    cond = pd.read_csv(os.path.join(RAW, "CondMeas.csv"))

    info = basic.merge(meta[["Sample_ID", "Lat", "Lon", "SamplingDepth",
                            "Source"]], on="Sample_ID", how="left")
    ksat = (cond[cond.MeasType == "KSAT"].groupby("Sample_ID").k.median()
            / 24.0)
    by_sample = {k: v for k, v in ret.groupby("Sample_ID")}

    rows = []
    skip = dict(points=0, texture=0, theta=0, fit=0, rmse=0)
    for r in info.itertuples():
        g = by_sample.get(r.Sample_ID)
        if g is None:
            skip["points"] += 1
            continue
        pf, theta = thin(g.pF.to_numpy(float), g.theta.to_numpy(float))
        if len(pf) < MIN_POINTS:
            skip["points"] += 1
            continue
        if theta.max() > MAX_THETA or theta.min() <= 0:
            skip["theta"] += 1
            continue
        h = 10.0 ** pf * CM_TO_KPA

        sand, silt, clay = r.Sand_USDA, r.Silt_USDA, r.Clay_USDA
        tot = sum(x for x in (sand, silt, clay) if np.isfinite(x))
        if not np.isfinite(tot) or not (95 <= tot <= 105):
            skip["texture"] += 1
            continue
        sand, silt, clay = (100 * sand / tot, 100 * silt / tot,
                            100 * clay / tot)

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

        oc = r.Corg if np.isfinite(r.Corg) else np.nan
        rows.append(dict(
            layer_id=f"HB_{r.Sample_ID}",
            # No profile field is published; samples sharing a site and depth
            # are replicates, so group on site to keep them out of each
            # other's leave-one-out folds.
            profile_id=f"HB_{r.Source}_{r.Lat:.4f}_{r.Lon:.4f}",
            texture_class=st.usda_class(sand, silt, clay),
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=sand, silt=silt, clay=clay,
            ksat_cmh=ksat.get(r.Sample_ID, np.nan),
            depth_cm=r.SamplingDepth, lat=r.Lat, lon=r.Lon,
            source_db="Hohenbrink", oc=oc, porosity=r.Porosity,
            rmse=rmse, n_points=len(pf), **fm))

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"samples read:    {len(info)}")
    print(f"reference layers: {len(out)}   (skipped: {skip})")
    print(f"  with ksat:      {out.ksat_cmh.notna().sum()}")
    print(f"  with oc:        {out.oc.notna().sum()}")
    print(f"  median points/curve after thinning: {int(out.n_points.median())}")
    print(f"  median fit RMSE: {out.rmse.median():.4f}")
    print(f"  profiles:        {out.profile_id.nunique()}")
    agree = (out.texture_class.to_numpy()
             == info.set_index("Sample_ID").TexClass_USDA.map(TEX)
             .reindex([i[3:] for i in out.layer_id]).to_numpy())
    print(f"  class from fractions matches published TexClass_USDA: "
          f"{agree.sum()}/{len(out)} ({agree.mean()*100:.1f} %)")
    print(out.texture_class.value_counts().to_string())


if __name__ == "__main__":
    main()
