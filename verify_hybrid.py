"""End-to-end check of the hybrid model as the tool actually runs it.

verify_gbm.py scored the classifier in isolation, predicting straight from
each target's fitted parameters. The shipped hybrid does something slightly
different: it averages predict_proba over the Monte Carlo draws of the vG fit
covariance, the same draws the kNN votes over. That is the right thing to do
-- it propagates fit uncertainty -- but it is not obviously free, so this
script measures the whole path through estimate().

Three things are checked:
  1. the hybrid reproduces the GBM's accuracy advantage over kNN once MC
     averaging is in the loop,
  2. Ksat, fractions and intervals are untouched (they must still come from
     the kNN),
  3. how often the kNN "second opinion" agrees, and whether agreement is
     actually informative about correctness.

Usage:  python verify_hybrid.py [n_per_class] [n_mc]
"""

import sys

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

import ksat_metrics as km
import swcc_texture as st
from verify_by_class import mcnemar
from verify_groups import GROUP
from verify_tau import prf
from verify_variants import COLS, ORDER

H = np.arange(0.0, 150.0 + 0.1, 5.0)
N_FOLDS = 5


def main():
    n_per_class = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    n_mc = int(sys.argv[2]) if len(sys.argv) > 2 else 40

    g = pd.read_csv(st.REFERENCE_CSV).reset_index(drop=True)
    g["profile_id"] = g.profile_id.fillna("solo_" + g.layer_id.astype(str))
    tsel = pd.concat([x.sample(min(len(x), n_per_class), random_state=0)
                      for _, x in g.groupby("texture_class")])
    tpos = np.array([g.index.get_loc(i) for i in tsel.index])
    truth = tsel.texture_class.to_numpy()
    obs = tsel.ksat_cmh.to_numpy(float)
    print(f"reference {len(g)} layers; targets {len(tsel)}; n_mc={n_mc}\n")

    n = len(tsel)
    p_knn, p_hyb, p_2nd = (np.empty(n, object) for _ in range(3))
    ks = {m: np.full(n, np.nan) for m in ("knn", "hyb")}
    lo = {m: np.full(n, np.nan) for m in ("knn", "hyb")}
    hi = {m: np.full(n, np.nan) for m in ("knn", "hyb")}

    for tr_idx, te_idx in GroupKFold(n_splits=N_FOLDS).split(
            g, groups=g.profile_id):
        sel = np.isin(tpos, te_idx)
        if not sel.any():
            continue
        train = g.iloc[tr_idx]
        ref = st.GshpReference(df=g[COLS].reset_index(drop=True))
        mask = np.ones(len(g), dtype=bool)
        mask[tr_idx] = False
        ref.set_excluded(mask)
        clf = st.TextureGBM(df=train)

        idx = np.where(sel)[0]
        for j, row in zip(idx, tsel[sel].itertuples()):
            theta = st.vg_theta(H, row.thetar, row.thetas, row.alpha_kpa, row.n)
            a = st.estimate(H, theta, ref=ref, n_mc=n_mc)
            b = st.estimate(H, theta, ref=ref, n_mc=n_mc, clf=clf)
            p_knn[j] = a["texture_class"]
            p_hyb[j] = b["texture_class"]
            p_2nd[j] = next(iter(b["knn_class_probabilities"]))
            for m, r in (("knn", a), ("hyb", b)):
                ks[m][j] = r["ksat"]["median_cmh"]
                lo[m][j] = r["ksat"]["p5_cmh"]
                hi[m][j] = r["ksat"]["p95_cmh"]

    p_knn, p_hyb, p_2nd = (x.astype(str) for x in (p_knn, p_hyb, p_2nd))
    print(f"{'model':<12s} {'overall':>8s} {'macroR':>7s} {'macroP':>7s} "
          f"{'macroF1':>8s} {'group':>7s}")
    print("-" * 52)
    for name, p in (("kNN", p_knn), ("hybrid", p_hyb)):
        R, P, F = prf(truth, p)
        grp = np.mean([GROUP[a] == GROUP[t] for a, t in zip(p, truth)])
        print(f"{name:<12s} {(p == truth).mean()*100:7.1f}% {R*100:6.1f}% "
              f"{P*100:6.1f}% {F*100:7.1f}% {grp*100:6.1f}%")
    a, b = p_knn == truth, p_hyb == truth
    print(f"\nhybrid vs kNN: {(b.mean()-a.mean())*100:+.1f}pp   "
          f"McNemar p = {mcnemar(a, b):.4g}")

    # 2. Ksat and intervals must be identical -- the GBM must not touch them.
    same = np.allclose(ks["knn"], ks["hyb"], equal_nan=True)
    print(f"\nKsat identical between kNN and hybrid: {same}  "
          f"(the GBM supplies only the class)")
    km.report("KSAT", [("kNN", km.ksat_scores(obs, ks["knn"], lo["knn"], hi["knn"])),
                       ("hybrid", km.ksat_scores(obs, ks["hyb"], lo["hyb"], hi["hyb"]))])

    # 3. Is the kNN second opinion informative?
    agree = p_hyb == p_2nd
    print(f"\nkNN second opinion agrees with the GBM on "
          f"{agree.mean()*100:.1f}% of soils")
    print(f"  hybrid accuracy when they agree:    {b[agree].mean()*100:5.1f}% "
          f"(n={agree.sum()})")
    print(f"  hybrid accuracy when they disagree: {b[~agree].mean()*100:5.1f}% "
          f"(n={(~agree).sum()})")

    print("\nPER-CLASS RECALL")
    print(f"{'class':<18s} {'n':>5s} {'kNN':>7s} {'hybrid':>7s} {'delta':>8s}")
    for c in ORDER:
        m = truth == c
        if m.any():
            print(f"{c:<18s} {m.sum():5d} {(p_knn[m]==c).mean()*100:6.1f}% "
                  f"{(p_hyb[m]==c).mean()*100:6.1f}% "
                  f"{((p_hyb[m]==c).mean()-(p_knn[m]==c).mean())*100:+7.1f}pp")


if __name__ == "__main__":
    main()
