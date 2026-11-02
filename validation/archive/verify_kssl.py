"""Does merging NCSS/KSSL into the reference improve the tool?

Paired profile-level leave-one-out over real GSHP soils, scoring BOTH
objectives:

  * texture  -- exact USDA class and aggregated group, compared with an exact
                McNemar test (the same soils are predicted under each
                reference, so the comparison is paired)
  * Ksat     -- log10 bias/RMSE, factor-of-N hit rates and coverage of the
                reported p5-p95 band (ksat_metrics)

Ksat matters especially here. KSSL distributes NO saturated conductivity at
all (see prepare_kssl.py), so every merged row is Ksat-blank. Only neighbours
carrying a measured Ksat contribute to the Ks estimate, so adding 2.5 k blank
rows can crowd Ksat-bearing neighbours out of the k=30 neighbourhood and make
Ks *worse* even if texture gets better. That is the trade-off this script
exists to measure.

Usage:  python verify_kssl.py [n_per_class] [n_mc] [--depth]
"""

import sys

import numpy as np
import pandas as pd

import ksat_metrics as km
import swcc_texture as st
from verify_by_class import mcnemar, predict_full
from verify_groups import GROUP, GROUP_ORDER
from verify_variants import COLS, ORDER

KSSL_CSV = "data/kssl_reference.csv"


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    use_depth = "--depth" in sys.argv
    n_per_class = int(args[0]) if args else 150
    n_mc = int(args[1]) if len(args) > 1 else 60

    gshp = pd.read_csv(st.REFERENCE_CSV)
    kssl = pd.read_csv(KSSL_CSV)
    merged = pd.concat([gshp[COLS], kssl[COLS]], ignore_index=True)

    pool = gshp[gshp.depth_cm.notna()]
    targets = pd.concat([g.sample(min(len(g), n_per_class), random_state=0)
                         for _, g in pool.groupby("texture_class")])
    truth = targets.texture_class.to_numpy()
    tgroup = np.array([GROUP[t] for t in truth])
    obs = targets.ksat_cmh.to_numpy(float)

    print(f"targets: {len(targets)} GSHP soils (stratified, <= {n_per_class}"
          f"/class); depth covariate: {use_depth}; n_mc={n_mc}")
    print(f"reference GSHP {len(gshp)} rows ({gshp.ksat_cmh.notna().sum()} "
          f"with Ksat)")
    print(f"reference GSHP+KSSL {len(merged)} rows "
          f"({merged.ksat_cmh.notna().sum()} with Ksat -- KSSL contributes 0, "
          f"so Ksat support falls from "
          f"{gshp.ksat_cmh.notna().mean()*100:.0f}% to "
          f"{merged.ksat_cmh.notna().mean()*100:.0f}% of rows)")
    print(f"KSSL adds {len(kssl)} layers\n")

    res = {}
    for name, ref in (("GSHP", gshp[COLS]), ("GSHP+KSSL", merged)):
        res[name] = predict_full(ref, targets, use_depth, n_mc)

    ok = {k: res[k][0] == truth for k in res}
    gok = {k: np.array([GROUP[p] == GROUP[t] for p, t in zip(res[k][0], truth)])
           for k in res}

    # ---- texture, by class -------------------------------------------------
    print("EXACT USDA CLASS")
    print(f"{'class':<18s} {'n':>5s} {'GSHP':>7s} {'+KSSL':>7s} {'delta':>8s} "
          f"{'p':>7s}")
    print("-" * 58)
    for cls in ORDER:
        m = truth == cls
        if not m.any():
            continue
        a, b = ok["GSHP"][m], ok["GSHP+KSSL"][m]
        print(f"{cls:<18s} {m.sum():5d} {a.mean()*100:6.1f}% {b.mean()*100:6.1f}% "
              f"{(b.mean()-a.mean())*100:+7.1f}pp {mcnemar(a, b):7.3f}")
    a, b = ok["GSHP"], ok["GSHP+KSSL"]
    print("-" * 58)
    print(f"{'ALL':<18s} {len(truth):5d} {a.mean()*100:6.1f}% {b.mean()*100:6.1f}% "
          f"{(b.mean()-a.mean())*100:+7.1f}pp {mcnemar(a, b):7.3f}")

    # ---- texture, by group -------------------------------------------------
    print("\nAGGREGATED GROUP")
    print(f"{'group':<18s} {'n':>5s} {'GSHP':>7s} {'+KSSL':>7s} {'delta':>8s} "
          f"{'p':>7s}")
    print("-" * 58)
    for g_ in GROUP_ORDER:
        m = tgroup == g_
        if not m.any():
            continue
        a, b = gok["GSHP"][m], gok["GSHP+KSSL"][m]
        print(f"{g_:<18s} {m.sum():5d} {a.mean()*100:6.1f}% {b.mean()*100:6.1f}% "
              f"{(b.mean()-a.mean())*100:+7.1f}pp {mcnemar(a, b):7.3f}")
    a, b = gok["GSHP"], gok["GSHP+KSSL"]
    print("-" * 58)
    print(f"{'ALL':<18s} {len(truth):5d} {a.mean()*100:6.1f}% {b.mean()*100:6.1f}% "
          f"{(b.mean()-a.mean())*100:+7.1f}pp {mcnemar(a, b):7.3f}")

    # ---- Ksat --------------------------------------------------------------
    km.report("KSAT (second objective)",
              [(k, km.ksat_scores(obs, res[k][1], res[k][2], res[k][3]))
               for k in ("GSHP", "GSHP+KSSL")])


if __name__ == "__main__":
    main()
