"""Tune the class-prior exponent tau.

Neighbour votes are weighted by (1/class_count)^tau. The tool has always used
tau = 1, a full uniform prior over the 12 USDA classes. That maximises
macro-recall but wrecks precision for rare classes: GSHP holds 32 silt layers
against 3,930 sand, so at tau = 1 one silt neighbour outvotes 123 sand
neighbours and silt is predicted about 2.7x more often than it occurs.

Recall alone cannot see this, which is why the selection criterion here is
macro-F1. Ksat is reported too, because tau reweights the same neighbours that
produce the Ks estimate.

CAVEAT ON PRECISION: precision depends on the class mix of the evaluation set.
Targets are stratified (as in every other script here) so that per-class recall
is comparable, but that makes the evaluation population roughly uniform, which
is the prior tau = 1 assumes. A natural-frequency evaluation would penalise
large tau harder. The stratified choice is therefore conservative -- it
flatters tau = 1.

Usage:  python verify_tau.py [n_per_class] [n_mc]
"""

import sys

import numpy as np
import pandas as pd

import ksat_metrics as km
import swcc_texture as st
from verify_groups import GROUP
from verify_variants import COLS, ORDER

TAUS = (0.0, 0.25, 0.5, 0.75, 1.0)
H = np.arange(0.0, 150.0 + 0.1, 5.0)


def prf(truth, preds):
    """Macro recall / precision / F1 over the classes present in truth."""
    recs, precs, f1s = [], [], []
    for c in ORDER:
        m, pm = truth == c, preds == c
        if not m.any():
            continue
        tp = int((m & pm).sum())
        r = tp / m.sum()
        p = tp / pm.sum() if pm.sum() else 0.0
        recs.append(r); precs.append(p)
        f1s.append(2 * p * r / (p + r) if (p + r) else 0.0)
    return np.mean(recs), np.mean(precs), np.mean(f1s)


def main():
    n_per_class = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    n_mc = int(sys.argv[2]) if len(sys.argv) > 2 else 60

    gshp = pd.read_csv(st.REFERENCE_CSV)
    targets = pd.concat([g.sample(min(len(g), n_per_class), random_state=0)
                         for _, g in gshp.groupby("texture_class")])
    truth = targets.texture_class.to_numpy()
    tgroup = np.array([GROUP[t] for t in truth])
    obs = targets.ksat_cmh.to_numpy(float)
    print(f"targets: {len(targets)} soils (stratified, <= {n_per_class}/class); "
          f"profile-level LOO; n_mc={n_mc}\n")

    print(f"{'tau':>5s} {'overall':>8s} {'macroR':>7s} {'macroP':>7s} "
          f"{'macroF1':>8s} {'group':>7s} {'silt P':>7s} {'silt pred':>10s}")
    print("-" * 64)
    ks_rows, best = [], None
    for tau in TAUS:
        ref = st.GshpReference(df=gshp[COLS].reset_index(drop=True), tau=tau)
        cls, med, lo, hi = [], [], [], []
        for row in targets.itertuples():
            theta = st.vg_theta(H, row.thetar, row.thetas, row.alpha_kpa, row.n)
            ref.set_excluded_profile(row.profile_id)
            r = st.estimate(H, theta, ref=ref, n_mc=n_mc)
            cls.append(r["texture_class"])
            k = r["ksat"]
            med.append(k["median_cmh"]); lo.append(k["p5_cmh"])
            hi.append(k["p95_cmh"])
        preds = np.array(cls)
        R, P, F = prf(truth, preds)
        overall = (preds == truth).mean()
        grp = np.mean([GROUP[p] == GROUP[t] for p, t in zip(preds, truth)])
        n_silt_pred = int((preds == "silt").sum())
        silt_p = (int(((truth == "silt") & (preds == "silt")).sum())
                  / n_silt_pred if n_silt_pred else float("nan"))
        print(f"{tau:5.2f} {overall*100:7.1f}% {R*100:6.1f}% {P*100:6.1f}% "
              f"{F*100:7.1f}% {grp*100:6.1f}% {silt_p*100:6.1f}% "
              f"{n_silt_pred:10d}")
        ks_rows.append((f"tau={tau:g}", km.ksat_scores(
            obs, np.array(med, float), np.array(lo, float),
            np.array(hi, float))))
        if best is None or F > best[1]:
            best = (tau, F)

    print(f"\n  true silt count in targets: {int((truth == 'silt').sum())}")
    print(f"  best macro-F1 at tau = {best[0]:g}")
    km.report("KSAT by tau", ks_rows)


if __name__ == "__main__":
    main()
