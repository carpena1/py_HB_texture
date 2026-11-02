"""Does the hybrid pass the absolute-skill tests that verify_skill.py applies?

verify_skill.py asked whether each class and group is acceptable in its own
right, under two null hypotheses:

  1. H0: recall equals the chance level (1/12 per class; group size / 12 per
     group, so Sandy 2/12, Loamy 4/12, Silty 3/12, Clayey 3/12).
     H1: recall greater. One-sided exact binomial.
  2. H0: recall <= 0.5, with --threshold 0.5.

This runs both tests for the kNN and the hybrid on IDENTICAL profile-grouped
folds, so the two are directly comparable.

NOTE ON COMPARABILITY WITH THE EARLIER verify_skill.py RUN: that one used
profile-level leave-one-out, which hides only the target's own profile. Here
both models must see the same training data, and the GBM has to be refitted
per fold, so the split is profile-grouped 5-fold -- each model sees 80 % of
the reference instead of ~100 %. kNN numbers are therefore a little lower
here than in the earlier run. The kNN-vs-hybrid comparison is unaffected,
because both are scored on exactly the same folds.

Usage:  python verify_skill_hybrid.py [n_per_class] [n_mc] [--threshold T]
"""

import sys

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

import swcc_texture as st
from verify_by_class import mcnemar
from verify_groups import GROUP, GROUP_ORDER
from verify_skill import assess, report
from verify_variants import COLS, ORDER

H = np.arange(0.0, 150.0 + 0.1, 5.0)
N_CLASSES = 12
N_FOLDS = 5


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
    n_mc = int(args[1]) if len(args) > 1 else 40

    g = pd.read_csv(st.REFERENCE_CSV).reset_index(drop=True)
    g["profile_id"] = g.profile_id.fillna("solo_" + g.layer_id.astype(str))
    tsel = pd.concat([x.sample(min(len(x), n_per_class), random_state=0)
                      for _, x in g.groupby("texture_class")])
    tpos = np.array([g.index.get_loc(i) for i in tsel.index])
    truth = tsel.texture_class.to_numpy()
    tgroup = np.array([GROUP[t] for t in truth])

    alpha = 0.05 / N_CLASSES
    print(f"targets: {len(tsel)} soils (stratified, <= {n_per_class}/class); "
          f"profile-grouped {N_FOLDS}-fold; n_mc={n_mc}")
    print(f"significance: alpha = 0.05/{N_CLASSES} = {alpha:.4f} "
          f"(Bonferroni over the 12 class tests)")
    if threshold is not None:
        print(f"absolute bar: H0 recall <= {threshold:.0%}")

    preds = {"kNN": np.empty(len(tsel), object),
             "hybrid": np.empty(len(tsel), object)}
    for tr_idx, te_idx in GroupKFold(n_splits=N_FOLDS).split(
            g, groups=g.profile_id):
        sel = np.isin(tpos, te_idx)
        if not sel.any():
            continue
        ref = st.GshpReference(df=g[COLS].reset_index(drop=True))
        mask = np.ones(len(g), dtype=bool)
        mask[tr_idx] = False
        ref.set_excluded(mask)
        clf = st.TextureGBM(df=g.iloc[tr_idx])
        for j, row in zip(np.where(sel)[0], tsel[sel].itertuples()):
            theta = st.vg_theta(H, row.thetar, row.thetas, row.alpha_kpa, row.n)
            preds["kNN"][j] = st.estimate(
                H, theta, ref=ref, n_mc=n_mc)["texture_class"]
            preds["hybrid"][j] = st.estimate(
                H, theta, ref=ref, n_mc=n_mc, clf=clf)["texture_class"]
    preds = {k: v.astype(str) for k, v in preds.items()}

    for name, p in preds.items():
        ok = p == truth
        gok = np.array([GROUP[a] == GROUP[t] for a, t in zip(p, truth)])
        print(f"\n{'=' * 78}\n{name}\n{'=' * 78}")

        rows = [assess(c, ok[truth == c].sum(), int((truth == c).sum()),
                       1 / N_CLASSES, threshold)
                for c in ORDER if (truth == c).any()]
        report("BY USDA CLASS (recall vs a uniform 1-of-12 guesser)", rows,
               threshold, alpha)

        print("\n  precision and F1 (same run):")
        print(f"  {'class':<18s} {'recall':>7s} {'precis.':>8s} {'F1':>7s} "
              f"{'predicted n':>12s}")
        for c in ORDER:
            if not (truth == c).any():
                continue
            r, pr, f1, npred = prec_f1(truth, p, c)
            print(f"  {c:<18s} {r*100:6.1f}% {pr*100:7.1f}% {f1*100:6.1f}% "
                  f"{npred:12d}")

        grows = [assess(gr, gok[tgroup == gr].sum(), int((tgroup == gr).sum()),
                        sum(1 for c in ORDER if GROUP[c] == gr) / N_CLASSES,
                        threshold)
                 for gr in GROUP_ORDER if (tgroup == gr).any()]
        report("BY AGGREGATED GROUP (chance = group size / 12)", grows,
               threshold, 0.05 / len(GROUP_ORDER),
               note="  NB: chance differs per group -- use the skill column.")

        o = assess("ALL", ok.sum(), len(ok), 1 / N_CLASSES, threshold)
        print(f"\n  overall exact-class recall {o['recall']*100:.1f}% "
              f"[{o['lo']*100:.1f},{o['hi']*100:.1f}], "
              f"skill {o['skill']*100:.1f}%, p(>chance) = {o['p_chance']:.2e}")

    # ---- side-by-side deltas ----------------------------------------------
    a, b = preds["kNN"], preds["hybrid"]
    print(f"\n{'=' * 78}\nkNN vs HYBRID, paired\n{'=' * 78}")
    print(f"{'class':<18s} {'n':>5s} {'kNN':>7s} {'hybrid':>7s} {'delta':>8s} "
          f"{'kNN F1':>7s} {'hyb F1':>7s} {'p':>7s}")
    print("-" * 72)
    for c in ORDER:
        m = truth == c
        if not m.any():
            continue
        _, _, f1a, _ = prec_f1(truth, a, c)
        _, _, f1b, _ = prec_f1(truth, b, c)
        ra, rb = (a[m] == c).mean(), (b[m] == c).mean()
        print(f"{c:<18s} {m.sum():5d} {ra*100:6.1f}% {rb*100:6.1f}% "
              f"{(rb-ra)*100:+7.1f}pp {f1a*100:6.1f}% {f1b*100:6.1f}% "
              f"{mcnemar(a[m]==c, b[m]==c):7.3f}")
    print("-" * 72)
    print(f"{'ALL':<18s} {len(truth):5d} {(a==truth).mean()*100:6.1f}% "
          f"{(b==truth).mean()*100:6.1f}% "
          f"{((b==truth).mean()-(a==truth).mean())*100:+7.1f}pp "
          f"{'':>15s} {mcnemar(a==truth, b==truth):7.3f}")

    print(f"\n{'group':<18s} {'n':>5s} {'kNN':>7s} {'hybrid':>7s} {'delta':>8s} "
          f"{'p':>7s}")
    print("-" * 58)
    ga = np.array([GROUP[x] == GROUP[t] for x, t in zip(a, truth)])
    gb = np.array([GROUP[x] == GROUP[t] for x, t in zip(b, truth)])
    for gr in GROUP_ORDER:
        m = tgroup == gr
        if m.any():
            print(f"{gr:<18s} {m.sum():5d} {ga[m].mean()*100:6.1f}% "
                  f"{gb[m].mean()*100:6.1f}% "
                  f"{(gb[m].mean()-ga[m].mean())*100:+7.1f}pp "
                  f"{mcnemar(ga[m], gb[m]):7.3f}")


if __name__ == "__main__":
    main()
