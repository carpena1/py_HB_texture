"""Is ~39 % a limit of the *method*, or of the *information*?

Nearest-neighbour lookup, two reference merges and three metric changes have
all landed in the same place. This script asks whether a strong discriminative
model does better on identical inputs and identical splits. If gradient
boosting cannot beat kNN, the ceiling is in the information that van Genuchten
parameters carry about texture, not in the lookup.

Compared on the same folds:
  * kNN          -- the tool, reference restricted to the training fold
  * kNN +qw      -- neighbours down-weighted by GSHP's own fit uncertainty
  * GBM          -- HistGradientBoostingClassifier on the same four features

Two split schemes:
  * profile-grouped 5-fold  -- comparable to the profile-level LOO used
                               everywhere else
  * source-blocked          -- leave one contributing laboratory out, the
                               honest out-of-lab estimate

The GBM's mean top probability is also reported. For a *calibrated* classifier
that would estimate the Bayes accuracy of this feature set, but the calibration
table printed below shows this model is overconfident (46.5 % mean confidence
against 39.4 % actual accuracy), so it must NOT be read as a ceiling. Use
verify_ceiling.py for a model-free bound.

Ksat is reported too: GBM regresses log10 Ksat, kNN uses its neighbour median.

Usage:  python verify_gbm.py [n_per_class] [n_mc] [--depth]
"""

import sys
import warnings

import numpy as np
import pandas as pd
from sklearn.ensemble import (HistGradientBoostingClassifier,
                              HistGradientBoostingRegressor)
from sklearn.model_selection import GroupKFold

import ksat_metrics as km
import swcc_texture as st
from verify_by_class import mcnemar
from verify_groups import GROUP
from verify_tau import prf
from verify_variants import ORDER

H = np.arange(0.0, 150.0 + 0.1, 5.0)
N_FOLDS = 5


def features(df, use_depth):
    cols = [np.log10(df.alpha_kpa), np.log10(df.n - 1.0), df.thetar, df.thetas]
    if use_depth:
        cols.append(np.log10(1.0 + df.depth_cm.clip(lower=0).fillna(0)))
    return np.column_stack(cols)


def knn_fold(gshp, train_idx, targets, n_mc, quality_weight):
    """kNN with the whole test fold hidden from the reference."""
    ref = st.GshpReference(df=gshp.reset_index(drop=True),
                           quality_weight=quality_weight)
    mask = np.ones(len(gshp), dtype=bool)
    mask[train_idx] = False           # hide everything not in train
    ref.set_excluded(mask)
    cls, med, lo, hi = [], [], [], []
    for row in targets.itertuples():
        theta = st.vg_theta(H, row.thetar, row.thetas, row.alpha_kpa, row.n)
        r = st.estimate(H, theta, ref=ref, n_mc=n_mc)
        cls.append(r["texture_class"])
        k = r["ksat"]
        med.append(k["median_cmh"]); lo.append(k["p5_cmh"]); hi.append(k["p95_cmh"])
    return np.array(cls), np.array(med, float), np.array(lo, float), np.array(hi, float)


def gbm_fold(gshp, train_idx, targets, use_depth):
    tr = gshp.iloc[train_idx]
    clf = HistGradientBoostingClassifier(
        max_iter=400, learning_rate=0.06, max_leaf_nodes=31,
        l2_regularization=1.0, early_stopping=True, validation_fraction=0.15,
        random_state=0, class_weight="balanced")
    clf.fit(features(tr, use_depth), tr.texture_class.to_numpy())
    Xt = features(targets, use_depth)
    proba = clf.predict_proba(Xt)
    preds = clf.classes_[proba.argmax(1)]

    # Ksat: regress log10 Ksat on the same features, train rows that have it.
    k = tr[tr.ksat_cmh.notna() & (tr.ksat_cmh > 0)]
    ks = np.full(len(targets), np.nan)
    if len(k) > 50:
        reg = HistGradientBoostingRegressor(
            max_iter=400, learning_rate=0.06, early_stopping=True,
            validation_fraction=0.15, random_state=0)
        reg.fit(features(k, use_depth), np.log10(k.ksat_cmh.to_numpy()))
        ks = 10.0 ** reg.predict(Xt)
    return preds, proba, clf.classes_, ks


def summarize(name, truth, preds, rows):
    R, P, F = prf(truth, preds)
    grp = np.mean([GROUP[p] == GROUP[t] for p, t in zip(preds, truth)])
    rows.append((name, (preds == truth).mean(), R, P, F, grp))


def table(title, rows):
    print(f"\n{title}")
    print(f"{'model':<16s} {'overall':>8s} {'macroR':>7s} {'macroP':>7s} "
          f"{'macroF1':>8s} {'group':>7s}")
    print("-" * 58)
    for n, o, R, P, F, g in rows:
        print(f"{n:<16s} {o*100:7.1f}% {R*100:6.1f}% {P*100:6.1f}% "
              f"{F*100:7.1f}% {g*100:6.1f}%")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    use_depth = "--depth" in sys.argv
    n_per_class = int(args[0]) if args else 150
    n_mc = int(args[1]) if len(args) > 1 else 40

    gshp = pd.read_csv(st.REFERENCE_CSV).reset_index(drop=True)
    if use_depth:
        gshp = gshp[gshp.depth_cm.notna()].reset_index(drop=True)
    # 168 layers carry no profile_id; treat each as its own profile so grouped
    # splitting keeps them separate rather than lumping them into one group.
    gshp["profile_id"] = gshp.profile_id.fillna(
        "solo_" + gshp.layer_id.astype(str))
    tsel = pd.concat([g.sample(min(len(g), n_per_class), random_state=0)
                      for _, g in gshp.groupby("texture_class")])
    tpos = np.array([gshp.index.get_loc(i) for i in tsel.index])
    truth = tsel.texture_class.to_numpy()
    obs = tsel.ksat_cmh.to_numpy(float)
    print(f"reference {len(gshp)} layers; targets {len(tsel)} "
          f"(stratified <= {n_per_class}/class); depth={use_depth}; n_mc={n_mc}")

    # ---------------- profile-grouped 5-fold --------------------------------
    gkf = GroupKFold(n_splits=N_FOLDS)
    pk, pkq, pg = (np.empty(len(tsel), object) for _ in range(3))
    pks = {k: np.full(len(tsel), np.nan) for k in ("knn", "knnq", "gbm")}
    plo = {k: np.full(len(tsel), np.nan) for k in ("knn", "knnq")}
    phi = {k: np.full(len(tsel), np.nan) for k in ("knn", "knnq")}
    proba_all = np.zeros((len(tsel), len(ORDER)))

    for tr_idx, te_idx in gkf.split(gshp, groups=gshp.profile_id):
        sel = np.isin(tpos, te_idx)          # targets inside this test fold
        if not sel.any():
            continue
        tg = tsel[sel]
        c, m, l, h = knn_fold(gshp, tr_idx, tg, n_mc, False)
        pk[sel], pks["knn"][sel], plo["knn"][sel], phi["knn"][sel] = c, m, l, h
        c, m, l, h = knn_fold(gshp, tr_idx, tg, n_mc, True)
        pkq[sel], pks["knnq"][sel], plo["knnq"][sel], phi["knnq"][sel] = c, m, l, h
        gp, proba, classes, gks = gbm_fold(gshp, tr_idx, tg, use_depth)
        pg[sel] = gp
        pks["gbm"][sel] = gks
        for j, cn in enumerate(classes):
            proba_all[np.where(sel)[0], ORDER.index(cn)] = proba[:, j]

    rows = []
    summarize("kNN", truth, pk.astype(str), rows)
    summarize("kNN +qual.wt", truth, pkq.astype(str), rows)
    summarize("GBM", truth, pg.astype(str), rows)
    table("PROFILE-GROUPED 5-FOLD", rows)
    base = pk.astype(str) == truth
    for nm, pp in (("kNN +qual.wt", pkq), ("GBM", pg)):
        b = pp.astype(str) == truth
        print(f"  {nm:<14s} vs kNN: {(b.mean()-base.mean())*100:+5.1f}pp  "
              f"McNemar p = {mcnemar(base, b):.4g}")

    # ---------------- per-class ---------------------------------------------
    print("\nPER-CLASS RECALL (profile-grouped)")
    # NB: the last column is the GBM's own mean top probability, NOT a Bayes
    # ceiling -- the calibration table below shows the model is overconfident.
    # For a model-free ceiling see verify_ceiling.py.
    print(f"{'class':<18s} {'n':>5s} {'kNN':>7s} {'GBM':>7s} {'delta':>8s} "
          f"{'GBM conf':>9s}")
    print("-" * 58)
    for c in ORDER:
        m = truth == c
        if not m.any():
            continue
        a = (pk[m].astype(str) == c).mean()
        b = (pg[m].astype(str) == c).mean()
        ceil = proba_all[m].max(1).mean()
        print(f"{c:<18s} {m.sum():5d} {a*100:6.1f}% {b*100:6.1f}% "
              f"{(b-a)*100:+7.1f}pp {ceil*100:7.1f}%")

    # ---------------- calibration + ceiling ---------------------------------
    conf = proba_all.max(1)
    hit = (pg.astype(str) == truth)
    print("\nGBM CALIBRATION (is the ceiling estimate trustworthy?)")
    print(f"{'confidence bin':<18s} {'n':>6s} {'mean conf':>10s} {'accuracy':>9s}")
    for lo_, hi_ in ((0, .2), (.2, .3), (.3, .4), (.4, .6), (.6, 1.01)):
        b = (conf >= lo_) & (conf < hi_)
        if b.sum():
            print(f"  [{lo_:.2f},{hi_:.2f})      {b.sum():6d} "
                  f"{conf[b].mean()*100:9.1f}% {hit[b].mean()*100:8.1f}%")
    print(f"\n  mean top probability  = {conf.mean()*100:.1f}%   "
          f"(estimated Bayes accuracy for these features)")
    print(f"  actual GBM accuracy   = {hit.mean()*100:.1f}%")

    # ---------------- source-blocked ----------------------------------------
    print("\n" + "=" * 58)
    rows = []
    pk2, pg2 = np.empty(len(tsel), object), np.empty(len(tsel), object)
    for s in gshp.source_db.unique():
        te = np.where(gshp.source_db.to_numpy() == s)[0]
        tr_idx = np.where(gshp.source_db.to_numpy() != s)[0]
        sel = np.isin(tpos, te)
        if not sel.any():
            continue
        tg = tsel[sel]
        pk2[sel] = knn_fold(gshp, tr_idx, tg, n_mc, False)[0]
        pg2[sel] = gbm_fold(gshp, tr_idx, tg, use_depth)[0]
    ok = np.array([x is not None for x in pk2])
    summarize("kNN", truth[ok], pk2[ok].astype(str), rows)
    summarize("GBM", truth[ok], pg2[ok].astype(str), rows)
    table(f"SOURCE-BLOCKED (n={ok.sum()})", rows)
    a2 = pk2[ok].astype(str) == truth[ok]
    b2 = pg2[ok].astype(str) == truth[ok]
    print(f"  {'GBM':<14s} vs kNN: {(b2.mean()-a2.mean())*100:+5.1f}pp  "
          f"McNemar p = {mcnemar(a2, b2):.4g}")

    # ---------------- Ksat ---------------------------------------------------
    km.report("KSAT (profile-grouped)", [
        ("kNN", km.ksat_scores(obs, pks["knn"], plo["knn"], phi["knn"])),
        ("kNN +qual.wt", km.ksat_scores(obs, pks["knnq"], plo["knnq"], phi["knnq"])),
        ("GBM regressor", km.ksat_scores(obs, pks["gbm"],
                                         np.full(len(obs), np.nan),
                                         np.full(len(obs), np.nan))),
    ])


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    main()
