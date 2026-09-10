"""What does the Hohenbrink et al. (2023) European dataset add?

Two questions, both reported for texture (exact class and group) and for Ksat:

 1. Does adding 560 German layers help EUROPEAN soils in the existing
    databases, and does it leave non-European soils alone? Targets are drawn
    from GSHP+KSSL only, so the gain cannot come from Hohenbrink predicting
    itself. The European/non-European split uses verify_region's bounding box.

 2. How well are the Hohenbrink soils themselves predicted from GSHP+KSSL,
    i.e. from laboratories that never measured them? This is the
    source-blocked view, and the honest estimate of what a user with a new
    German core should expect.

Both models are scored. Profile-grouped folds are used rather than
leave-one-out so that the GBM's training frame and the kNN's reference are
always the same, and neither can see the target.

Usage:  python verify_hohenbrink.py [n_per_class] [n_mc]
"""

import sys

import numpy as np
import pandas as pd

import ksat_metrics as km
import swcc_texture as st
from verify_by_class import mcnemar
from verify_groups import GROUP, GROUP_ORDER
from verify_region import europe_mask
from verify_tau import prf
from verify_variants import ORDER

H = np.concatenate([[0.0], np.logspace(-1, np.log10(1500.0), 25)])
N_FOLDS = 5


def curve(row):
    return st.vg_theta(H, row.thetar, row.thetas, row.alpha_kpa, row.n)


def predict(train_df, targets, n_mc, hybrid, folds=None):
    """If folds is None the whole train_df is the reference and targets are
    assumed disjoint from it; otherwise profile-grouped folds are removed
    from both the kNN reference and the GBM training frame."""
    cls = np.empty(len(targets), dtype=object)
    med = np.full(len(targets), np.nan)
    lo = np.full(len(targets), np.nan)
    hi = np.full(len(targets), np.nan)

    def do(train, sel):
        ref = st.GshpReference(df=train.reset_index(drop=True))
        clf = st.TextureGBM(df=train.reset_index(drop=True)) if hybrid else None
        for j in sel:
            row = targets.iloc[j]
            r = st.estimate(H, curve(row), ref=ref, n_mc=n_mc, clf=clf)
            cls[j] = r["texture_class"]
            k = r["ksat"]
            med[j], lo[j], hi[j] = k["median_cmh"], k["p5_cmh"], k["p95_cmh"]

    if folds is None:
        do(train_df, range(len(targets)))
    else:
        tp = targets.profile_id.to_numpy()
        for held in folds:
            do(train_df[~train_df.profile_id.isin(held)],
               np.where(np.isin(tp, list(held)))[0])
    return cls, med, lo, hi


def score(name, truth, pred, obs, med, lo, hi):
    R, P, F = prf(truth, pred)
    grp = np.mean([GROUP[a] == GROUP[t] for a, t in zip(pred, truth)])
    print(f"{name:<26s} {(pred == truth).mean()*100:7.1f}% {R*100:6.1f}% "
          f"{F*100:7.1f}% {grp*100:6.1f}%")
    return km.ksat_scores(obs, med, lo, hi)


def block(title, truth, obs, res):
    print(f"\n=== {title} (n={len(truth)}) ===")
    print(f"{'variant':<26s} {'exact':>8s} {'macroR':>7s} {'macroF1':>8s} "
          f"{'group':>7s}")
    ks = []
    for name, (p, m, l, h) in res.items():
        ks.append((name, score(name, truth, p, obs, m, l, h)))
    keys = list(res)
    for i in range(0, len(keys) - 1, 2):
        a, b = keys[i], keys[i + 1]
        ca, cb = res[a][0] == truth, res[b][0] == truth
        ga = np.array([GROUP[x] == GROUP[t] for x, t in zip(res[a][0], truth)])
        gb = np.array([GROUP[x] == GROUP[t] for x, t in zip(res[b][0], truth)])
        print(f"  {b} vs {a}: class {(cb.mean()-ca.mean())*100:+.1f}pp "
              f"p={mcnemar(ca, cb):.4f}   group "
              f"{(gb.mean()-ga.mean())*100:+.1f}pp p={mcnemar(ga, gb):.4f}")
    km.report(f"KSAT -- {title}", ks)

    print(f"\n  per-class recall ({title})")
    hdr = "".join(f"{k:>16s}" for k in res)
    print(f"  {'class':<18s} {'n':>5s}{hdr}")
    for cl in ORDER:
        m = truth == cl
        if not m.any():
            continue
        vals = "".join(f"{(res[k][0][m] == cl).mean()*100:15.1f}%" for k in res)
        print(f"  {cl:<18s} {m.sum():5d}{vals}")
    print(f"\n  per-group accuracy ({title})")
    print(f"  {'group':<18s} {'n':>5s}{hdr}")
    for g in GROUP_ORDER:
        m = np.array([GROUP[t] for t in truth]) == g
        if not m.any():
            continue
        vals = "".join(
            f"{np.mean([GROUP[x] == GROUP[t] for x, t in zip(res[k][0][m], truth[m])])*100:15.1f}%"
            for k in res)
        print(f"  {g:<18s} {m.sum():5d}{vals}")


def main():
    n_per_class = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    n_mc = int(sys.argv[2]) if len(sys.argv) > 2 else 40

    base = pd.concat([st.load_reference_df("gshp"), st.load_reference_df("kssl")],
                     ignore_index=True)
    hb = st.load_reference_df("hohenbrink")
    both = pd.concat([base, hb], ignore_index=True)
    for d in (base, hb, both):
        d["profile_id"] = d.profile_id.fillna("solo_" + d.layer_id.astype(str))

    eu = europe_mask(base)
    print(f"reference without Hohenbrink: {len(base)}  with: {len(both)}")
    print(f"European layers in base: {int(eu.sum())}  "
          f"(+{len(hb)} from Hohenbrink)")
    print(f"targets stratified {n_per_class}/class; profile-grouped "
          f"{N_FOLDS}-fold; n_mc={n_mc}\n")

    rng = np.random.default_rng(0)
    for label, mask in [("EUROPEAN targets", eu), ("NON-EUROPEAN targets", ~eu)]:
        pool = base[mask]
        tg = pd.concat([g.sample(min(len(g), n_per_class), random_state=0)
                        for _, g in pool.groupby("texture_class")])
        truth = tg.texture_class.to_numpy()
        obs = tg.ksat_cmh.to_numpy(float)
        profs = np.array(sorted(tg.profile_id.unique()))
        rng.shuffle(profs)
        folds = [set(x) for x in np.array_split(profs, N_FOLDS)]
        res = {}
        for tag, hybrid in (("kNN", False), ("hybrid", True)):
            res[f"{tag} base"] = predict(base, tg, n_mc, hybrid, folds)
            res[f"{tag} +Hohenbrink"] = predict(both, tg, n_mc, hybrid, folds)
        order = ["kNN base", "kNN +Hohenbrink", "hybrid base",
                 "hybrid +Hohenbrink"]
        block(label, truth, obs, {k: res[k] for k in order})

    # --- Hohenbrink itself, with and without its own source in the reference
    # The paired arms answer a question the source-blocked figure alone
    # cannot: is the Ksat bias a property of these soils, or simply of a
    # reference whose Ksat population sits 1.5 orders of magnitude lower?
    tg = hb
    truth = tg.texture_class.to_numpy()
    obs = tg.ksat_cmh.to_numpy(float)
    hprofs = np.array(sorted(tg.profile_id.unique()))
    np.random.default_rng(1).shuffle(hprofs)
    hfolds = [set(x) for x in np.array_split(hprofs, N_FOLDS)]
    block("HOHENBRINK soils: reference without vs with their own source",
          truth, obs,
          {"kNN, source-blocked": predict(base, tg, n_mc, False),
           "kNN, +own source": predict(both, tg, n_mc, False, hfolds),
           "hybrid, source-blocked": predict(base, tg, n_mc, True),
           "hybrid, +own source": predict(both, tg, n_mc, True, hfolds)})


if __name__ == "__main__":
    main()
