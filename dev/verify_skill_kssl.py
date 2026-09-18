"""Detailed per-class and per-group skill for GSHP alone vs GSHP+KSSL (kNN).

verify_kssl.py reported only the headline delta (+0.5 pp, p=0.523). This gives
the full breakdown the skill tests use: recall with a Wilson interval,
chance-corrected skill, precision/F1, and both null hypotheses -- recall vs
chance, and recall vs an absolute 50 % bar.

Design notes:
  * targets are GSHP soils only; KSSL is purely extra reference material and
    is never a target,
  * profile-level leave-one-out (the target's own profile is hidden), which is
    the tool's real operating condition and matches verify_kssl.py. KSSL
    profile ids are prefixed "KSSL_" so they never collide with GSHP ones and
    are never excluded,
  * kNN only -- no GBM -- since this is about what the reference table
    contributes.

Because this is leave-one-out rather than the 5-fold used in
verify_skill_hybrid.py, the GSHP-only column here sits roughly 1.5 points
higher than the kNN column there. Both are correct; they answer different
questions.

Usage:  python verify_skill_kssl.py [n_per_class] [n_mc] [--threshold T]
"""

import sys

import numpy as np
import pandas as pd

import ksat_metrics as km
import swcc_texture as st
from verify_by_class import mcnemar, predict_full
from verify_groups import GROUP, GROUP_ORDER
from verify_skill import assess, report
from verify_variants import COLS, ORDER

N_CLASSES = 12
KSSL_CSV = "data/kssl_reference.csv"


def prec_f1(truth, preds, cls):
    m, pm = truth == cls, preds == cls
    tp = int((m & pm).sum())
    r = tp / m.sum() if m.any() else float("nan")
    p = tp / pm.sum() if pm.sum() else float("nan")
    f = (2 * p * r / (p + r) if pm.sum() and m.any() and (p + r) > 0
         else float("nan"))
    return r, p, f, int(pm.sum())


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    threshold = None
    for i, a in enumerate(sys.argv[1:], start=1):
        if a.startswith("--threshold"):
            threshold = (float(a.split("=")[1]) if "=" in a
                         else float(sys.argv[i + 1]))
    n_per_class = int(args[0]) if args else 150
    n_mc = int(args[1]) if len(args) > 1 else 50

    gshp = pd.read_csv(st.REFERENCE_CSV)
    kssl = pd.read_csv(KSSL_CSV)
    merged = pd.concat([gshp[COLS], kssl[COLS]], ignore_index=True)

    targets = pd.concat([g.sample(min(len(g), n_per_class), random_state=0)
                         for _, g in gshp.groupby("texture_class")])
    truth = targets.texture_class.to_numpy()
    tgroup = np.array([GROUP[t] for t in truth])
    obs = targets.ksat_cmh.to_numpy(float)

    alpha = 0.05 / N_CLASSES
    print(f"targets: {len(targets)} GSHP soils (stratified <= {n_per_class}"
          f"/class); profile-level LOO; kNN; n_mc={n_mc}")
    print(f"GSHP {len(gshp)} layers; GSHP+KSSL {len(merged)} "
          f"(KSSL adds {len(kssl)}, all with the 0.95*porosity wet-end anchor)")
    print(f"significance: alpha = 0.05/{N_CLASSES} = {alpha:.4f} (Bonferroni)")
    if threshold is not None:
        print(f"absolute bar: H0 recall <= {threshold:.0%}")

    res = {}
    for name, ref in (("GSHP", gshp[COLS]), ("GSHP+KSSL", merged)):
        res[name] = predict_full(ref, targets, False, n_mc)

    for name in res:
        p = res[name][0]
        ok = p == truth
        gok = np.array([GROUP[a] == GROUP[t] for a, t in zip(p, truth)])
        print(f"\n{'=' * 78}\n{name}\n{'=' * 78}")

        rows = [assess(c, ok[truth == c].sum(), int((truth == c).sum()),
                       1 / N_CLASSES, threshold)
                for c in ORDER if (truth == c).any()]
        report("BY USDA CLASS", rows, threshold, alpha)

        print("\n  precision and F1 (same run):")
        print(f"  {'class':<18s} {'recall':>7s} {'precis.':>8s} {'F1':>7s} "
              f"{'predicted n':>12s}")
        for c in ORDER:
            if (truth == c).any():
                r, pr, f1, npred = prec_f1(truth, p, c)
                print(f"  {c:<18s} {r*100:6.1f}% {pr*100:7.1f}% {f1*100:6.1f}% "
                      f"{npred:12d}")

        grows = [assess(gr, gok[tgroup == gr].sum(), int((tgroup == gr).sum()),
                        sum(1 for c in ORDER if GROUP[c] == gr) / N_CLASSES,
                        threshold)
                 for gr in GROUP_ORDER if (tgroup == gr).any()]
        report("BY AGGREGATED GROUP (chance = group size / 12)", grows,
               threshold, 0.05 / len(GROUP_ORDER))

        o = assess("ALL", ok.sum(), len(ok), 1 / N_CLASSES, threshold)
        print(f"\n  overall {o['recall']*100:.1f}% "
              f"[{o['lo']*100:.1f},{o['hi']*100:.1f}], "
              f"skill {o['skill']*100:.1f}%")

    # ---- paired deltas -----------------------------------------------------
    a, b = res["GSHP"][0], res["GSHP+KSSL"][0]
    print(f"\n{'=' * 78}\nPAIRED: GSHP vs GSHP+KSSL\n{'=' * 78}")
    print(f"{'class':<18s} {'n':>5s} {'GSHP':>7s} {'+KSSL':>7s} {'delta':>8s} "
          f"{'F1 GSHP':>8s} {'F1 +K':>7s} {'p':>7s}")
    print("-" * 74)
    for c in ORDER:
        m = truth == c
        if not m.any():
            continue
        _, _, f1a, _ = prec_f1(truth, a, c)
        _, _, f1b, _ = prec_f1(truth, b, c)
        print(f"{c:<18s} {m.sum():5d} {(a[m]==c).mean()*100:6.1f}% "
              f"{(b[m]==c).mean()*100:6.1f}% "
              f"{((b[m]==c).mean()-(a[m]==c).mean())*100:+7.1f}pp "
              f"{f1a*100:7.1f}% {f1b*100:6.1f}% "
              f"{mcnemar(a[m]==c, b[m]==c):7.3f}")
    print("-" * 74)
    print(f"{'ALL':<18s} {len(truth):5d} {(a==truth).mean()*100:6.1f}% "
          f"{(b==truth).mean()*100:6.1f}% "
          f"{((b==truth).mean()-(a==truth).mean())*100:+7.1f}pp "
          f"{'':>16s} {mcnemar(a==truth, b==truth):7.3f}")

    ga = np.array([GROUP[x] == GROUP[t] for x, t in zip(a, truth)])
    gb = np.array([GROUP[x] == GROUP[t] for x, t in zip(b, truth)])
    print(f"\n{'group':<18s} {'n':>5s} {'GSHP':>7s} {'+KSSL':>7s} {'delta':>8s} "
          f"{'p':>7s}")
    print("-" * 58)
    for gr in GROUP_ORDER:
        m = tgroup == gr
        if m.any():
            print(f"{gr:<18s} {m.sum():5d} {ga[m].mean()*100:6.1f}% "
                  f"{gb[m].mean()*100:6.1f}% "
                  f"{(gb[m].mean()-ga[m].mean())*100:+7.1f}pp "
                  f"{mcnemar(ga[m], gb[m]):7.3f}")
    print("-" * 58)
    print(f"{'ALL':<18s} {len(truth):5d} {ga.mean()*100:6.1f}% "
          f"{gb.mean()*100:6.1f}% {(gb.mean()-ga.mean())*100:+7.1f}pp "
          f"{mcnemar(ga, gb):7.3f}")

    km.report("KSAT", [(n, km.ksat_scores(obs, res[n][1], res[n][2], res[n][3]))
                       for n in res])


if __name__ == "__main__":
    main()
