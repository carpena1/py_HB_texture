"""IvG fits (ivg.py) for every reference layer, cached for verify_ivg.py.

    python ivg_fits.py [n_procs]

Layers whose source publishes measured retention points are refitted to
those points, over their whole suction range (Arizona's run to ~2e5 kPa,
beyond the 1500 kPa its reference row was fitted to). The rest -- GSHP's
databases, UNSODA, and the Yellow River set, whose points are rebuilt from
the authors' fits -- publish only vG parameters, so their curve is sampled
at h = 0 and 12 log-spaced suctions from 0.1 to 1500 kPa and refitted: the
IvG parameters then describe the same vG curve and carry no dry-end
information of their own. `points` records which.

Three fits per layer: the paper's thetar bound (thetas/20, columns ivg_*),
the looser absolute bound of the Arizona spreadsheet (0.03, ivgA_*), and the
paper's bound with m free of n (ivg_*_m). Each row also carries theta at
fixed heads (THETA_HEADS_CM) from the first two and from the reference's own
vG parameters, for the fixed-head predictors. Writes figures/cache/ivg_fits.csv (git-ignored: it includes
EU-HYDI layers).
"""

import importlib
import os
import sys
import warnings
from multiprocessing import Pool

import numpy as np
import pandas as pd

import ivg
import swcc_texture as st

OUT = os.path.join("figures", "cache", "ivg_fits.csv")
CM_TO_KPA = 0.0980665
# theta_s (h = 0) and the middle and dry heads the brief lists, in cm
THETA_HEADS_CM = np.array([0.0, 50, 100, 330, 1000, 5000, 15000, 100000])
SYNTHETIC_KPA = np.r_[0.0, np.logspace(-1, np.log10(1500.0), 12)]
MEASURED = {"kssl": {}, "hohenbrink": {}, "euhydi": {},
            "babaeian_zanjanrood": {}, "babaeian_az": {"max_kpa": np.inf},
            "armas": {}, "willard": {}, "boorowa": {}}


def points():
    pts = {}
    for name, kw in MEASURED.items():
        try:
            pts.update(importlib.import_module("prepare_" + name)
                       .measured_points(**kw))
        except (FileNotFoundError, OSError, SystemExit):
            print(f"  {name}: no measured points here", file=sys.stderr)
    return pts


def fit_one(args):
    lid, h, theta, kind = args
    warnings.filterwarnings("ignore")
    out = dict(layer_id=lid, points=kind, n_pts=len(h),
               h_max_kpa=float(np.max(h)))
    for stem, suffix, bound, free in (("ivg", "", "thetas20", False),
                                      ("ivgA", "", "abs", False),
                                      ("ivg", "_m", "thetas20", True)):
        try:
            tr, ts, a, n, m, rmse = ivg.fit_ivg(h, theta, free_m=free,
                                                bound=bound)
        except Exception:
            tr = ts = a = n = m = rmse = np.nan
        out.update({f"{stem}_thetar{suffix}": tr, f"{stem}_thetas{suffix}": ts,
                    f"{stem}_alpha_kpa{suffix}": a, f"{stem}_n{suffix}": n,
                    f"{stem}_rmse{suffix}": rmse})
        if free:
            if not np.isfinite(rmse):
                # too few points for five parameters: keep the tied-m fit,
                # as free_m.py does
                tr, ts, a, n = (out[f"{stem}_thetar"], out[f"{stem}_thetas"],
                                out[f"{stem}_alpha_kpa"], out[f"{stem}_n"])
                m, rmse = 1.0 - 1.0 / n, out[f"{stem}_rmse"]
                out.update({f"{stem}_thetar{suffix}": tr,
                            f"{stem}_thetas{suffix}": ts,
                            f"{stem}_alpha_kpa{suffix}": a,
                            f"{stem}_n{suffix}": n, f"{stem}_rmse{suffix}": rmse})
            out[f"{stem}_m{suffix}"] = m
        else:
            th = ivg.ivg_theta(THETA_HEADS_CM * CM_TO_KPA, tr, ts, a, n)
            out.update({f"{stem}_th{int(c)}": v
                        for c, v in zip(THETA_HEADS_CM, th)})
    return out


def main():
    procs = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    ref = st.load_reference_df("all").drop_duplicates("layer_id")
    # the Boorowa test set is not in "all"; fit it too
    try:
        ref = pd.concat([ref, st.load_reference_df("boorowa")],
                        ignore_index=True)
    except Exception:
        pass
    pts = points()
    jobs = []
    for r in ref.itertuples():
        if r.layer_id in pts:
            h, th = pts[r.layer_id]
            jobs.append((r.layer_id, np.asarray(h, float),
                         np.asarray(th, float), "measured"))
        else:
            th = st.vg_theta(SYNTHETIC_KPA, r.thetar, r.thetas, r.alpha_kpa,
                             r.n)
            jobs.append((r.layer_id, SYNTHETIC_KPA, th, "synthetic"))
    print(f"{len(jobs)} layers: "
          f"{sum(j[3] == 'measured' for j in jobs)} measured, "
          f"{sum(j[3] == 'synthetic' for j in jobs)} synthetic", flush=True)
    with Pool(procs) as p:
        rows = p.map(fit_one, jobs, chunksize=50)
    out = pd.DataFrame(rows).set_index("layer_id")
    r = ref.set_index("layer_id")
    vg = st.vg_theta(THETA_HEADS_CM[None, :] * CM_TO_KPA,
                     r.thetar.to_numpy()[:, None], r.thetas.to_numpy()[:, None],
                     r.alpha_kpa.to_numpy()[:, None], r.n.to_numpy()[:, None])
    for j, c in enumerate(THETA_HEADS_CM):
        out.loc[r.index, f"vg_th{int(c)}"] = vg[:, j]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    out.to_csv(OUT)
    print(f"wrote {OUT}; fit failures: "
          f"{int(out.ivg_rmse.isna().sum())} (thetas/20), "
          f"{int(out.ivgA_rmse.isna().sum())} (0.03), "
          f"{int(out.ivg_rmse_m.isna().sum())} (free m, after fallback)")
    print(out.groupby("points")[["ivg_rmse", "ivgA_rmse",
                                 "ivg_rmse_m"]].median())


if __name__ == "__main__":
    main()
