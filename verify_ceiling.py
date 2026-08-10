"""A model-free estimate of the accuracy ceiling.

verify_gbm.py reports the GBM's mean top probability, but that model is
overconfident, so it cannot be read as a Bayes-accuracy estimate. This script
uses the classical Cover & Hart (1967) result instead, which needs no model
and no calibration: asymptotically the 1-nearest-neighbour error rate brackets
the Bayes error,

    err_1NN / 2  <=  err_Bayes  <=  err_1NN

so Bayes accuracy lies between 1 - err_1NN and 1 - err_1NN/2. Whatever
classifier we build on these four van Genuchten features, it cannot beat the
upper end of that bracket.

Splits are profile-grouped, so a layer is never its own neighbour and never
matched against a sibling horizon from the same site. A source-blocked variant
is also reported, since that is the situation a real user is in.

THREE CAVEATS, in order of importance:

1. Only the OVERALL bracket is rigorous. Cover & Hart bounds the total error
   rate; applying it to a single class's recall is a heuristic, which is why a
   few per-class results land outside their own bracket. Read the per-class
   columns as an indicator of relative separability, not as a bound.

2. The bounds are asymptotic. With finite data the measured 1NN error exceeds
   its asymptotic value, so "Bayes accuracy >= 1 - err_1NN" stays valid but
   the upper end is an estimate rather than a guarantee -- the true ceiling
   could sit somewhat above it.

3. Bayes error depends on the class prior. GSHP is 39 % sand while the tool
   imposes a uniform prior, so both are reported. Note the balanced reference
   is also smaller (capped per class), which confounds prior with reference
   size; compare the two with that in mind.

Usage:  python verify_ceiling.py [n_per_class] [cap_per_class]
"""

import sys

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

import swcc_texture as st
from verify_groups import GROUP
from verify_variants import ORDER

N_FOLDS = 5


def feats(df):
    x = np.column_stack([np.log10(df.alpha_kpa), np.log10(df.n - 1.0),
                         df.thetar, df.thetas])
    return x


def one_nn(train, test, mean, std):
    a = (feats(train) - mean) / std
    b = (feats(test) - mean) / std
    out = np.empty(len(b), dtype=object)
    cls = train.texture_class.to_numpy()
    for i in range(len(b)):
        d = ((a - b[i]) ** 2).sum(1)
        out[i] = cls[d.argmin()]
    return out


def bayes_bracket(err_1nn, n_classes):
    """Cover & Hart bounds on Bayes accuracy from the 1NN error rate.

    The multi-class form is tighter than the familiar two-class one:

        err_1NN <= err_B * (2 - M/(M-1) * err_B)

    Inverting for the smallest err_B consistent with the observed err_1NN
    gives the upper bound on Bayes accuracy; err_B <= err_1NN gives the lower.
    """
    M = n_classes
    disc = 1.0 - (M / (M - 1.0)) * err_1nn
    e_min = ((M - 1.0) / M) * (1.0 - np.sqrt(disc)) if disc > 0 else \
        (M - 1.0) / M
    return 1.0 - err_1nn, 1.0 - e_min


def bracket(err, label, n_classes=12):
    lo, hi = bayes_bracket(err, n_classes)
    print(f"  {label:<28s} 1NN accuracy {(1-err)*100:5.1f}%  ->  "
          f"Bayes accuracy in [{lo*100:.1f}%, {hi*100:.1f}%]")
    return lo, hi


def run(gshp, targets, tpos, groups, label):
    mean, std = feats(gshp).mean(0), feats(gshp).std(0)
    truth = targets.texture_class.to_numpy()
    preds = np.empty(len(targets), dtype=object)
    gkf = GroupKFold(n_splits=N_FOLDS)
    for tr, te in gkf.split(gshp, groups=groups):
        sel = np.isin(tpos, te)
        if sel.any():
            preds[sel] = one_nn(gshp.iloc[tr], targets[sel], mean, std)
    ok = preds == truth
    print(f"\n{label}")
    bracket(1 - ok.mean(), "overall", 12)
    grp = np.mean([GROUP[p] == GROUP[t] for p, t in zip(preds, truth)])
    bracket(1 - grp, "aggregated group", 4)
    print(f"\n  {'class':<18s} {'n':>5s} {'1NN':>7s} "
          f"{'Bayes accuracy bracket':>24s}")
    for c in ORDER:
        m = truth == c
        if not m.any():
            continue
        e = 1 - (preds[m] == c).mean()
        lo, hi = bayes_bracket(e, 12)
        print(f"  {c:<18s} {m.sum():5d} {(1-e)*100:6.1f}%   "
              f"[{lo*100:5.1f}%, {hi*100:5.1f}%]")
    return ok


def main():
    n_per_class = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    cap = int(sys.argv[2]) if len(sys.argv) > 2 else 300

    g = pd.read_csv(st.REFERENCE_CSV).reset_index(drop=True)
    g["profile_id"] = g.profile_id.fillna("solo_" + g.layer_id.astype(str))

    for name, ref in (("NATURAL PRIOR (GSHP as-is, 39 % sand)", g),
                      (f"BALANCED PRIOR (<= {cap}/class, the prior the tool "
                       f"assumes)",
                       pd.concat([x.sample(min(len(x), cap), random_state=0)
                                  for _, x in g.groupby("texture_class")])
                       .reset_index(drop=True))):
        t = pd.concat([x.sample(min(len(x), n_per_class), random_state=1)
                       for _, x in ref.groupby("texture_class")])
        tpos = np.array([ref.index.get_loc(i) for i in t.index])
        print("\n" + "=" * 66)
        print(name + f"   (reference {len(ref)} layers, targets {len(t)})")
        run(ref, t, tpos, ref.profile_id.to_numpy(),
            "profile-grouped 5-fold")
        run(ref, t, tpos, ref.source_db.to_numpy(),
            "source-blocked (grouped by contributing laboratory)")


if __name__ == "__main__":
    main()
