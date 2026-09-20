"""Neighbour-mean particle fractions against sand/clay regressors, on
outside sources.

The kNN reports fractions as the weighted mean of its neighbours'. The
alternative is to predict sand and clay with gradient-boosted regressors and
take silt = 100 - sand - clay (swcc_texture.FractionGBM), which keeps the sum
constraint and drops silt, the weakest of the three.

On targets drawn from the reference the regressors win (verify_alt_predictors
.py: with depth and bulk density, silt 11.5 vs 12.7 points, sand 12.6 vs
13.4). This script is the test that decided against them: every external set
predicted as a NEW SOURCE, its own table held out of the reference, from its
measured points, with both sets of fractions taken from the same run so the
comparison is paired.

The regressors lose here, and lose on clay in particular: boosting shrinks
extreme clays towards the mean, while the class-frequency weighting lets the
neighbour mean reach them. Laikipia (51 % clay on average) is the clearest
case.

Usage:  python verify_fractions.py [--cov] [n_mc]
        --cov  supply depth and bulk density, as the command line does
"""

import importlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import swcc_texture as st

SETS = {"babaeian_az": "Arizona", "willard": "Laikipia", "boorowa": "Boorowa",
        "babaeian_zanjanrood": "Zanjanrood", "armas": "Canary"}


def main():
    use_cov = "--cov" in sys.argv
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    n_mc = int(args[0]) if args else 60

    g = st.load_reference_df("merged")
    print(f"{'set':12s} {'n':>4s}  " + "".join(f"{c + ' knn->gbm':>19s}"
                                               for c in ("sand", "silt", "clay"))
          + "   mean")
    tot = {"knn": [], "gbm": []}
    for name, label in SETS.items():
        m = importlib.import_module("prepare_" + name)
        try:
            tab = pd.read_csv(m.OUT)
        except FileNotFoundError:
            print(f"{label:12s}    -   ({m.OUT} not built here)")
            continue
        ref_df = g[~g.layer_id.isin(tab.layer_id)].reset_index(drop=True)
        ref = st.GshpReference(df=ref_df, use_bd=use_cov)
        clf = st.TextureGBM(df=ref_df, fractions=True,
                            covariates=["depth_cm", "bd"] if use_cov else [])
        pts = m.measured_points()
        rows = {"knn": [], "gbm": []}
        for r in tab.itertuples():
            if r.layer_id not in pts:
                continue
            if use_cov and not (np.isfinite(r.bd) and np.isfinite(r.depth_cm)):
                continue
            h, th = pts[r.layer_id]
            res = st.estimate(h, th, ref=ref, clf=clf, n_mc=n_mc,
                              sample_type=r.sample_type,
                              depth=r.depth_cm if use_cov else None,
                              bulk_density=r.bd if use_cov else None)
            truth = np.array([r.sand, r.silt, r.clay])
            for k, key in (("knn", "fractions"), ("gbm", "fractions_gbm")):
                rows[k].append(np.abs(np.array(
                    [res[key][c] for c in ("sand", "silt", "clay")]) - truth))
        if not rows["knn"]:
            print(f"{label:12s}    0   (no soil carries both depth and bulk "
                  f"density)")
            continue
        a, b = np.array(rows["knn"]), np.array(rows["gbm"])
        tot["knn"].append(a)
        tot["gbm"].append(b)
        cells = "".join(f"  {a[:, j].mean():8.1f} ->{b[:, j].mean():6.1f}"
                        for j in range(3))
        print(f"{label:12s} {len(a):4d} {cells}   {a.mean():5.2f} -> "
              f"{b.mean():5.2f}")

    A, B = np.vstack(tot["knn"]), np.vstack(tot["gbm"])
    print(f"\nall {len(A)} soils{' with depth and bulk density' if use_cov else ''}: "
          f"mean |error| {A.mean():.2f} -> {B.mean():.2f} points "
          f"(Wilcoxon p={wilcoxon(A.mean(1), B.mean(1)).pvalue:.4f}); "
          f"regressors better for {np.mean(B.mean(1) < A.mean(1)) * 100:.0f} % "
          f"of soils")
    for j, c in enumerate(("sand", "silt", "clay")):
        print(f"  {c:5s} {A[:, j].mean():5.2f} -> {B[:, j].mean():5.2f}  "
              f"(p={wilcoxon(A[:, j], B[:, j]).pvalue:.4f})")


if __name__ == "__main__":
    main()
