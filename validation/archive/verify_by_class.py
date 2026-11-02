"""Accuracy broken down by USDA texture class and by aggregated texture group,
globally (all soils), with and without the depth covariate.

Reports BOTH objectives: texture classification and Ksat. Ksat is scored in
log10 space by ksat_metrics (bias, RMSE, factor-of-N hit rates, and coverage
of the reported p5-p95 band).

Targets are a stratified sample of real GSHP soils that all carry a measured
depth, so the with- and without-depth runs are scored on exactly the same
soils and the comparison is paired. Each target is predicted with its whole
profile removed from the reference (profile-level leave-one-out).

Usage:  python verify_by_class.py [n_per_class] [n_mc]
"""

import sys
from math import comb

import numpy as np
import pandas as pd

import ksat_metrics as km
import swcc_texture as st
from verify_groups import GROUP, GROUP_ORDER
from verify_variants import COLS, ORDER

H = np.arange(0.0, 150.0 + 0.1, 5.0)


def mcnemar(a_ok, b_ok):
    b = int(np.sum(a_ok & ~b_ok))
    c = int(np.sum(~a_ok & b_ok))
    n = b + c
    if n == 0:
        return 1.0
    lo = min(b, c)
    return min(1.0, sum(comb(n, i) for i in range(lo + 1)) / 2 ** n * 2)


def predict_full(ref_df, targets, use_depth, n_mc):
    """Return predicted class plus the Ksat estimate and its p5-p95 band.

    Ksat is a second objective alongside texture, and it comes out of the same
    neighbour lookup, so scoring it here costs nothing extra.
    """
    ref = st.GshpReference(df=ref_df.reset_index(drop=True),
                           use_depth=use_depth)
    cls, med, lo, hi = [], [], [], []
    for row in targets.itertuples():
        theta = st.vg_theta(H, row.thetar, row.thetas, row.alpha_kpa, row.n)
        ref.set_excluded_profile(row.profile_id)
        res = st.estimate(H, theta, ref=ref, n_mc=n_mc,
                          depth=row.depth_cm if use_depth else None)
        cls.append(res["texture_class"])
        k = res["ksat"]
        med.append(k["median_cmh"]); lo.append(k["p5_cmh"]); hi.append(k["p95_cmh"])
    return (np.array(cls), np.array(med, float), np.array(lo, float),
            np.array(hi, float))


def predict(ref_df, targets, use_depth, n_mc):
    return predict_full(ref_df, targets, use_depth, n_mc)[0]


def main():
    n_per_class = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    n_mc = int(sys.argv[2]) if len(sys.argv) > 2 else 60

    gshp = pd.read_csv(st.REFERENCE_CSV)
    pool = gshp[gshp.depth_cm.notna()]
    targets = pd.concat([g.sample(min(len(g), n_per_class), random_state=0)
                         for _, g in pool.groupby("texture_class")])
    truth = targets.texture_class.to_numpy()
    print(f"global targets: {len(targets)} soils (stratified, <= {n_per_class}"
          f"/class, all with measured depth); reference {len(gshp)}; "
          f"n_mc={n_mc}\n")

    full = {d: predict_full(gshp[COLS], targets, d, n_mc) for d in (False, True)}
    preds = {d: full[d][0] for d in (False, True)}
    ok = {d: preds[d] == truth for d in (False, True)}
    grp_ok = {d: np.array([GROUP[p] == GROUP[t]
                           for p, t in zip(preds[d], truth)])
              for d in (False, True)}

    # ---- by USDA class -----------------------------------------------------
    print("BY USDA TEXTURE CLASS (exact-class accuracy)")
    print(f"{'class':<18s} {'n':>5s} {'no depth':>9s} {'+depth':>8s} "
          f"{'delta':>8s} {'p':>7s}")
    print("-" * 60)
    for cls in ORDER:
        m = truth == cls
        if not m.any():
            continue
        a, b = ok[False][m], ok[True][m]
        print(f"{cls:<18s} {m.sum():5d} {a.mean()*100:8.1f}% "
              f"{b.mean()*100:7.1f}% {(b.mean()-a.mean())*100:+7.1f}pp "
              f"{mcnemar(a, b):7.3f}")
    a, b = ok[False], ok[True]
    print("-" * 60)
    print(f"{'ALL (macro avg)':<18s} {'':5s} "
          f"{np.mean([ok[False][truth==c].mean() for c in ORDER if (truth==c).any()])*100:7.1f}% "
          f"{np.mean([ok[True][truth==c].mean() for c in ORDER if (truth==c).any()])*100:7.1f}%")
    print(f"{'ALL (overall)':<18s} {len(truth):5d} {a.mean()*100:8.1f}% "
          f"{b.mean()*100:7.1f}% {(b.mean()-a.mean())*100:+7.1f}pp "
          f"{mcnemar(a, b):7.3f}")

    # ---- by aggregated group ----------------------------------------------
    print("\nBY AGGREGATED GROUP (correct group, i.e. near-misses forgiven)")
    print(f"{'group':<18s} {'n':>5s} {'no depth':>9s} {'+depth':>8s} "
          f"{'delta':>8s} {'p':>7s}")
    print("-" * 60)
    tg = np.array([GROUP[t] for t in truth])
    for g in GROUP_ORDER:
        m = tg == g
        if not m.any():
            continue
        a, b = grp_ok[False][m], grp_ok[True][m]
        print(f"{g:<18s} {m.sum():5d} {a.mean()*100:8.1f}% "
              f"{b.mean()*100:7.1f}% {(b.mean()-a.mean())*100:+7.1f}pp "
              f"{mcnemar(a, b):7.3f}")
    a, b = grp_ok[False], grp_ok[True]
    print("-" * 60)
    print(f"{'ALL (overall)':<18s} {len(truth):5d} {a.mean()*100:8.1f}% "
          f"{b.mean()*100:7.1f}% {(b.mean()-a.mean())*100:+7.1f}pp "
          f"{mcnemar(a, b):7.3f}")

    # ---- Ksat (second objective) ------------------------------------------
    obs = targets.ksat_cmh.to_numpy(float)
    rows = []
    for d in (False, True):
        _, med, lo, hi = full[d]
        rows.append(("+depth" if d else "no depth",
                     km.ksat_scores(obs, med, lo, hi)))
    km.report("KSAT, ALL TARGETS", rows)

    for d in (False, True):
        _, med, lo, hi = full[d]
        rows = []
        for cls in ORDER:
            m = truth == cls
            if m.any() and np.isfinite(obs[m]).any():
                rows.append((cls, km.ksat_scores(obs[m], med[m], lo[m], hi[m])))
        km.report(f"KSAT BY CLASS ({'+depth' if d else 'no depth'})", rows)


if __name__ == "__main__":
    main()
