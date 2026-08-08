"""Large-sample, paired leave-one-out comparison of the reference variants.

verify_variants.py scores each variant on 12 soils per benchmark, which cannot
separate differences of one or two soils from noise. This script runs the same
four variants over a large stratified sample of real soils, leaving each target
out of the reference, and compares the variants *pairwise on the same soils*
(a McNemar-style test), which is what makes small differences interpretable.

Usage:  python verify_variants_large.py [n_per_class] [n_mc]
"""

import sys

import numpy as np
import pandas as pd

import swcc_texture as st
from verify_groups import GROUP
from verify_variants import COLS, ORDER

H = np.arange(0.0, 150.0 + 0.1, 5.0)


def sample_targets(df, n_per_class, seed=0):
    """Stratified sample: up to n_per_class soils per USDA class, restricted
    to soils that have a depth (so every variant sees the same targets)."""
    d = df[df.depth_cm.notna()]
    parts = [g.sample(min(len(g), n_per_class), random_state=seed)
             for _, g in d.groupby("texture_class")]
    return pd.concat(parts).reset_index(drop=True)


def predict_all(ref_df, use_depth, targets, n_mc):
    ref = st.GshpReference(df=ref_df.reset_index(drop=True),
                           use_depth=use_depth)
    preds = []
    for row in targets.itertuples():
        theta = st.vg_theta(H, row.thetar, row.thetas, row.alpha_kpa, row.n)
        ref.set_excluded(row.layer_id)
        depth = row.depth_cm if use_depth else None
        preds.append(st.estimate(H, theta, ref=ref, n_mc=n_mc,
                                 depth=depth)["texture_class"])
    ref.set_excluded(None)
    return np.array(preds)


def mcnemar(a_ok, b_ok):
    """Exact binomial p-value for paired disagreements (two-sided)."""
    from math import comb
    b = int(np.sum(a_ok & ~b_ok))     # a right, b wrong
    c = int(np.sum(~a_ok & b_ok))     # b right, a wrong
    n = b + c
    if n == 0:
        return b, c, 1.0
    lo = min(b, c)
    p = sum(comb(n, i) for i in range(lo + 1)) / 2 ** n * 2
    return b, c, min(1.0, p)


def main():
    n_per_class = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    n_mc = int(sys.argv[2]) if len(sys.argv) > 2 else 50

    gshp = pd.read_csv(st.REFERENCE_CSV)[COLS]
    unsoda = pd.read_csv("data/unsoda_reference.csv")[COLS]
    merged = pd.concat([gshp, unsoda], ignore_index=True)

    tg = sample_targets(gshp, n_per_class)
    tu = sample_targets(unsoda, n_per_class)
    print(f"targets: {len(tg)} GSHP soils, {len(tu)} UNSODA soils "
          f"(stratified, <= {n_per_class}/class; n_mc={n_mc})")
    print(f"reference: GSHP {len(gshp)}, merged {len(merged)}\n")

    variants = [("GSHP", gshp, False), ("GSHP+depth", gshp, True),
                ("GSHP+UNSODA", merged, False),
                ("GSHP+UNSODA+depth", merged, True)]

    results = {}
    for vname, vdf, vdepth in variants:
        ok_cls, ok_grp = [], []
        for targets in (tg, tu):
            preds = predict_all(vdf, vdepth, targets, n_mc)
            truth = targets.texture_class.to_numpy()
            ok_cls.append(preds == truth)
            ok_grp.append(np.array([GROUP[p] == GROUP[t]
                                    for p, t in zip(preds, truth)]))
        results[vname] = (np.concatenate(ok_cls), np.concatenate(ok_grp),
                          [a.mean() for a in ok_cls])
        c, g, per = results[vname]
        print(f"{vname:<20s} class {c.mean()*100:5.1f} %   "
              f"group {g.mean()*100:5.1f} %   "
              f"(GSHP {per[0]*100:.0f} % / UNSODA {per[1]*100:.0f} % class)")

    print("\n=== paired comparisons vs the GSHP baseline (same soils) ===")
    base_c, base_g, _ = results["GSHP"]
    n = len(base_c)
    print(f"{'variant':<20s} {'level':<6s} {'delta':>8s} {'b':>4s} {'c':>4s} "
          f"{'p':>8s}")
    for vname in ["GSHP+depth", "GSHP+UNSODA", "GSHP+UNSODA+depth"]:
        vc, vg, _ = results[vname]
        for level, base, var in [("class", base_c, vc), ("group", base_g, vg)]:
            b, c, p = mcnemar(base, var)
            delta = (var.mean() - base.mean()) * 100
            print(f"{vname:<20s} {level:<6s} {delta:+7.1f}pp {b:4d} {c:4d} "
                  f"{p:8.3f}")
    print(f"\nn = {n} paired soils. b = baseline right/variant wrong, "
          f"c = variant right/baseline wrong.\np < 0.05 means the difference "
          f"is unlikely to be chance.")


if __name__ == "__main__":
    main()
