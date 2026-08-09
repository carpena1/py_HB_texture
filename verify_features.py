"""Does matching on the retention curve beat matching on the vG parameters?

The tool has always measured neighbour distance in standardized
(log10 alpha, log10(n-1), thetar, thetas). But alpha and n trade off along the
fit ridge, so that space distorts curve similarity: two soils with nearly
identical curves can sit far apart, and two soils close in parameter space can
have visibly different curves.

feature_mode="curve" instead matches on theta evaluated at eight fixed matric
potentials (swcc_texture.CURVE_HEADS: 1, 3, 10, 33, 100, 330, 1000, 1500 kPa),
which is close to an L2 distance between the curves themselves on a log-h grid
and includes the two classical anchors, field capacity and wilting point.

Crossed with tau, since both levers reweight the same neighbourhood and could
interact. Paired: every configuration predicts the same targets under
profile-level leave-one-out. Texture and Ksat are both reported.

Usage:  python verify_features.py [n_per_class] [n_mc]
"""

import sys

import numpy as np
import pandas as pd

import ksat_metrics as km
import swcc_texture as st
from verify_by_class import mcnemar
from verify_groups import GROUP
from verify_tau import prf
from verify_variants import COLS, ORDER

H = np.arange(0.0, 150.0 + 0.1, 5.0)
CONFIGS = [("vg", 1.0), ("curve", 1.0), ("curve_white", 1.0),
           ("vg", 0.75), ("curve_white", 0.75)]


def run(gshp, targets, mode, tau, n_mc):
    ref = st.GshpReference(df=gshp[COLS].reset_index(drop=True), tau=tau,
                           feature_mode=mode)
    cls, med, lo, hi = [], [], [], []
    for row in targets.itertuples():
        theta = st.vg_theta(H, row.thetar, row.thetas, row.alpha_kpa, row.n)
        ref.set_excluded_profile(row.profile_id)
        r = st.estimate(H, theta, ref=ref, n_mc=n_mc)
        cls.append(r["texture_class"])
        k = r["ksat"]
        med.append(k["median_cmh"]); lo.append(k["p5_cmh"]); hi.append(k["p95_cmh"])
    return (np.array(cls), np.array(med, float), np.array(lo, float),
            np.array(hi, float))


def main():
    n_per_class = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    n_mc = int(sys.argv[2]) if len(sys.argv) > 2 else 50

    gshp = pd.read_csv(st.REFERENCE_CSV)
    targets = pd.concat([g.sample(min(len(g), n_per_class), random_state=0)
                         for _, g in gshp.groupby("texture_class")])
    truth = targets.texture_class.to_numpy()
    obs = targets.ksat_cmh.to_numpy(float)
    print(f"targets: {len(targets)} soils (stratified, <= {n_per_class}/class); "
          f"profile-level LOO; n_mc={n_mc}")
    print(f"curve heads (kPa): {list(st.CURVE_HEADS)}\n")

    res, ks_rows = {}, []
    print(f"{'features':>9s} {'tau':>5s} {'overall':>8s} {'macroR':>7s} "
          f"{'macroP':>7s} {'macroF1':>8s} {'group':>7s}")
    print("-" * 54)
    for mode, tau in CONFIGS:
        out = run(gshp, targets, mode, tau, n_mc)
        res[(mode, tau)] = out
        preds = out[0]
        R, P, F = prf(truth, preds)
        grp = np.mean([GROUP[p] == GROUP[t] for p, t in zip(preds, truth)])
        print(f"{mode:>9s} {tau:5.2f} {(preds==truth).mean()*100:7.1f}% "
              f"{R*100:6.1f}% {P*100:6.1f}% {F*100:7.1f}% {grp*100:6.1f}%")
        ks_rows.append((f"{mode} tau={tau:g}",
                        km.ksat_scores(obs, out[1], out[2], out[3])))

    base = res[("vg", 1.0)][0] == truth
    print("\npaired McNemar against the current default (vg, tau=1):")
    for cfg in CONFIGS[1:]:
        b = res[cfg][0] == truth
        print(f"  {cfg[0]:>5s} tau={cfg[1]:<5g} {(b.mean()-base.mean())*100:+6.1f}pp"
              f"   p = {mcnemar(base, b):.4f}")

    print("\nper-class, current default vs best feature mode at the same tau:")
    alt = res[("curve_white", 1.0)][0]
    cur = res[("vg", 1.0)][0]
    print(f"{'class':<18s} {'n':>5s} {'vg':>7s} {'curve':>7s} {'delta':>8s}")
    for c in ORDER:
        m = truth == c
        if not m.any():
            continue
        a, b = (cur[m] == c).mean(), (alt[m] == c).mean()
        print(f"{c:<18s} {m.sum():5d} {a*100:6.1f}% {b*100:6.1f}% "
              f"{(b-a)*100:+7.1f}pp")

    km.report("KSAT", ks_rows)


if __name__ == "__main__":
    main()
