"""Is the classifier's performance ACCEPTABLE for each texture class and each
aggregated group -- judged on its own, separately for the with- and
without-depth configurations?

This asks a different question from verify_by_class.py. That script tests
whether adding depth *changes* accuracy (a paired McNemar test). Here each
configuration is judged in absolute terms.

Two null hypotheses are offered:

  1. Skill over chance (default).
     H0: recall for this class equals what a chance classifier would achieve.
     H1: recall is greater (one-sided).
     The chance classifier guesses uniformly over the 12 USDA classes, which
     matches the uniform prior the tool imposes through inverse class-frequency
     weighting. So chance recall is 1/12 = 8.3 % for every class. For an
     aggregated group the chance level is (number of classes in the group)/12,
     i.e. Sandy 2/12, Loamy 4/12, Silty 3/12, Clayey 3/12 -- group results are
     NOT comparable to each other without this correction.

  2. An absolute bar (--threshold T).
     H0: recall <= T.   H1: recall > T.
     Use this when "acceptable" means a fixed standard rather than "beats
     chance", e.g. --threshold 0.5.

Both are exact one-sided binomial tests. A Wilson 95 % confidence interval is
reported so the precision of each estimate is visible, and a skill score
(recall - chance)/(1 - chance) gives the chance-corrected performance.
Precision and F1 are also reported: high recall with poor precision is not an
acceptable classifier, and recall alone cannot show that.

Usage:  python verify_skill.py [n_per_class] [n_mc] [--threshold T]
"""

import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest

import swcc_texture as st
from verify_by_class import predict
from verify_groups import GROUP, GROUP_ORDER
from verify_variants import COLS, ORDER

N_CLASSES = 12


def wilson(k, n, z=1.96):
    """Wilson score interval for a binomial proportion."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z ** 2 / n
    c = (p + z ** 2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / d
    return max(0.0, c - h), min(1.0, c + h)


def assess(label, correct, n, chance, threshold):
    k = int(correct)
    rec = k / n if n else float("nan")
    lo, hi = wilson(k, n)
    p_chance = binomtest(k, n, chance, alternative="greater").pvalue
    skill = (rec - chance) / (1 - chance) if n else float("nan")
    out = dict(label=label, n=n, recall=rec, lo=lo, hi=hi, chance=chance,
               skill=skill, p_chance=p_chance)
    if threshold is not None:
        out["p_thresh"] = binomtest(k, n, threshold,
                                    alternative="greater").pvalue
    return out


def report(title, rows, threshold, alpha, note=""):
    print(f"\n{title}")
    if note:
        print(note)
    hdr = (f"{'':<18s} {'n':>4s} {'recall':>7s} {'95% CI':>15s} "
           f"{'chance':>7s} {'skill':>7s} {'p(>chance)':>11s}")
    if threshold is not None:
        hdr += f" {'p(>'+f'{threshold:.0%}'+')':>10s}"
    print(hdr)
    print("-" * (len(hdr) + 2))
    for r in rows:
        line = (f"{r['label']:<18s} {r['n']:4d} {r['recall']*100:6.1f}% "
                f"[{r['lo']*100:4.1f},{r['hi']*100:5.1f}] "
                f"{r['chance']*100:6.1f}% {r['skill']*100:6.1f}% "
                f"{r['p_chance']:11.2e}")
        if threshold is not None:
            line += f" {r['p_thresh']:10.3f}"
        flag = "" if r["p_chance"] < alpha else "   <-- NOT above chance"
        print(line + flag)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    threshold = None
    for a in sys.argv[1:]:
        if a.startswith("--threshold"):
            threshold = float(a.split("=")[1]) if "=" in a else \
                float(sys.argv[sys.argv.index(a) + 1])
    n_per_class = int(args[0]) if args else 150
    n_mc = int(args[1]) if len(args) > 1 else 60

    gshp = pd.read_csv(st.REFERENCE_CSV)
    pool = gshp[gshp.depth_cm.notna()]
    targets = pd.concat([g.sample(min(len(g), n_per_class), random_state=0)
                         for _, g in pool.groupby("texture_class")])
    truth = targets.texture_class.to_numpy()
    tgroup = np.array([GROUP[t] for t in truth])

    # Bonferroni over the 12 class tests reported per configuration.
    alpha = 0.05 / N_CLASSES
    print(f"targets: {len(targets)} soils (stratified, <= {n_per_class}/class, "
          f"all with measured depth); reference {len(gshp)}; n_mc={n_mc}")
    print(f"significance: alpha = 0.05/{N_CLASSES} = {alpha:.4f} "
          f"(Bonferroni over the 12 class tests)")

    for use_depth in (False, True):
        cfg = "WITH depth" if use_depth else "WITHOUT depth"
        preds = predict(gshp[COLS], targets, use_depth, n_mc)
        ok = preds == truth
        gok = np.array([GROUP[p] == GROUP[t] for p, t in zip(preds, truth)])

        print(f"\n{'=' * 78}\n{cfg}\n{'=' * 78}")

        rows = []
        for cls in ORDER:
            m = truth == cls
            if m.any():
                rows.append(assess(cls, ok[m].sum(), int(m.sum()),
                                   1 / N_CLASSES, threshold))
        report("BY USDA CLASS (recall vs a uniform 1-of-12 guesser)", rows,
               threshold, alpha)

        # Precision / F1 per class -- recall alone cannot show over-prediction.
        print("\n  precision and F1 (same run):")
        print(f"  {'class':<18s} {'recall':>7s} {'precis.':>8s} {'F1':>7s} "
              f"{'predicted n':>12s}")
        for cls in ORDER:
            m = truth == cls
            pm = preds == cls
            if not m.any():
                continue
            tp = int((m & pm).sum())
            rec = tp / m.sum()
            prec = tp / pm.sum() if pm.sum() else float("nan")
            f1 = (2 * prec * rec / (prec + rec)
                  if pm.sum() and (prec + rec) > 0 else float("nan"))
            print(f"  {cls:<18s} {rec*100:6.1f}% {prec*100:7.1f}% "
                  f"{f1*100:6.1f}% {int(pm.sum()):12d}")

        grows = []
        for g in GROUP_ORDER:
            m = tgroup == g
            n_in = sum(1 for c in ORDER if GROUP[c] == g)
            if m.any():
                grows.append(assess(g, gok[m].sum(), int(m.sum()),
                                    n_in / N_CLASSES, threshold))
        report("BY AGGREGATED GROUP (chance = group size / 12)", grows,
               threshold, 0.05 / len(GROUP_ORDER),
               note="  NB: chance differs per group, so raw recalls are not "
                    "comparable across groups -- use the skill column.")

        overall = assess("ALL", ok.sum(), len(ok), 1 / N_CLASSES, threshold)
        print(f"\n  overall exact-class recall {overall['recall']*100:.1f}% "
              f"[{overall['lo']*100:.1f},{overall['hi']*100:.1f}], "
              f"skill {overall['skill']*100:.1f}%, "
              f"p(>chance) = {overall['p_chance']:.2e}")


if __name__ == "__main__":
    main()
