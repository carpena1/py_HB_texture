"""Alternative predictors: IvG parameters, fixed-head water contents, and
predicting silt by difference.

The shipped tool matches soils on four Mualem-constrained van Genuchten
parameters. Three alternatives are tested here, all on identical targets and
folds, so only the predictors change:

  vg        log10 alpha, log10(n-1), thetar, thetas        (the shipped set)
  vg_m      the same with m freed from n                   (free_m.py)
  ivg       IvG (ivg.py) with the paper's thetar bound
  ivgA      IvG with the looser bound thetar <= 0.03
  ivg_m     IvG, paper's bound, m freed from n
  th_ivg    theta at 0, 50, 100, 330, 1000, 5000, 15000 and 100000 cm,
            from the IvG fit
  th_vg     the same heads from the vG fit, which separates what the fixed
            heads do from what the IvG model does

Each arm is run with the curve alone and with depth and bulk density added,
as the command line does with --depth and --bulk-density.

Texture class comes from the same gradient-boosted classifier the tool uses;
particle fractions come two ways -- the neighbour-weighted mean the tool
uses, and gradient-boosted regressors for sand and clay with silt taken as
100 - sand - clay (the third item in the brief) -- and Ks from the
neighbours' weighted median, as the tool does.

Unlike the shipped pipeline this runs one fit per soil rather than Monte
Carlo draws around it, for every arm alike: the comparison between arms is
like for like, but the absolute numbers are not directly comparable with
README's.

Needs figures/cache/ivg_fits.csv (python ivg_fits.py).

With --measured-only the reference and the targets keep just the layers
whose source publishes measured retention points, so no IvG fit in the run
is a refit of a published vG curve. That is the brief's option of leaving
the fitted-parameter-only databases out.

With the design "ceiling" the script reports, for each predictor set, the
Cover & Hart (1967) bracket on Bayes accuracy that verify_ceiling.py reports
for the shipped features: no classifier built on that set can do better than
the top of its bracket. That answers whether a predictor set carries more
texture information at all, separately from whether our models exploit it.

Usage:  python verify_ivg.py [profile|source|az|ceiling] [n_per_class]
                             [--measured-only]
"""

import os
import sys

import numpy as np
import pandas as pd

import swcc_texture as st
from verify_ceiling import bayes_bracket
from verify_common import GROUP, mcnemar

CACHE = os.path.join("figures", "cache", "ivg_fits.csv")
HEADS = [0, 50, 100, 330, 1000, 5000, 15000, 100000]
K = 30
N_FOLDS = 5
ARMS = {
    "vg": ["log_alpha", "log_n1", "thetar", "thetas"],
    "vg_m": ["log_alpha_m", "log_n1_m", "thetar_m", "thetas_m", "m"],
    "ivg": ["ivg_log_alpha", "ivg_log_n1", "ivg_thetar", "ivg_thetas"],
    "ivgA": ["ivgA_log_alpha", "ivgA_log_n1", "ivgA_thetar", "ivgA_thetas"],
    "ivg_m": ["ivg_log_alpha_m", "ivg_log_n1_m", "ivg_thetar_m",
              "ivg_thetas_m", "ivg_m_m"],
    "th_ivg": [f"ivg_th{h}" for h in HEADS],
    "th_vg": [f"vg_th{h}" for h in HEADS],
}
COVS = ["log_depth", "bd"]


def features(df):
    """Every arm's predictors as columns of one frame."""
    f = pd.DataFrame(index=df.index)
    for out, a, n, tr, ts in (
            ("", "alpha_kpa", "n", "thetar", "thetas"),
            ("_m", "alpha_kpa_m", "n_m", "thetar_m", "thetas_m"),
            ("ivg_", "ivg_alpha_kpa", "ivg_n", "ivg_thetar", "ivg_thetas"),
            ("ivgA_", "ivgA_alpha_kpa", "ivgA_n", "ivgA_thetar", "ivgA_thetas"),
            ("ivg_m", "ivg_alpha_kpa_m", "ivg_n_m", "ivg_thetar_m",
             "ivg_thetas_m")):
        # "ivg_m" names the free-m IvG arm's columns ivg_log_alpha_m, ...
        pre, suf = (out, "") if out.endswith("_") else ("", out)
        if out == "ivg_m":
            pre, suf = "ivg_", "_m"
        f[pre + "log_alpha" + suf] = np.log10(df[a])
        f[pre + "log_n1" + suf] = np.log10(df[n] - 1.0)
        f[pre + "thetar" + suf] = df[tr].to_numpy()
        f[pre + "thetas" + suf] = df[ts].to_numpy()
    f["m"] = df["m"].to_numpy()
    f["ivg_m_m"] = df["ivg_m_m"].to_numpy()
    for h in HEADS:
        f[f"ivg_th{h}"] = df[f"ivg_th{h}"].to_numpy()
        f[f"vg_th{h}"] = df[f"vg_th{h}"].to_numpy()
    f["log_depth"] = np.log10(1.0 + df["depth_cm"].clip(lower=0))
    return f          # bulk density is used as it stands, from the reference


def gbm(random_state=0, regress=False):
    from sklearn.ensemble import (HistGradientBoostingClassifier,
                                  HistGradientBoostingRegressor)
    kw = dict(max_iter=400, learning_rate=0.06, max_leaf_nodes=31,
              l2_regularization=1.0, early_stopping=True,
              validation_fraction=0.15, random_state=random_state)
    if regress:
        return HistGradientBoostingRegressor(**kw)
    return HistGradientBoostingClassifier(class_weight="balanced", **kw)


def knn(ztrain, ztest, train, k=K):
    """Neighbour fractions and Ks for each test row, as the tool weights
    them: 1/d^2 times the class-frequency weight."""
    d2 = ((ztest[:, None, :] - ztrain[None, :, :]) ** 2).sum(axis=2)
    idx = np.argpartition(d2, k, axis=1)[:, :k]
    frac = np.empty((len(ztest), 3))
    ks = np.full(len(ztest), np.nan)
    cw = train["class_weight"].to_numpy()
    fr = train[["sand", "silt", "clay"]].to_numpy()
    ksat = train["ksat_cmh"].to_numpy()
    for i, ix in enumerate(idx):
        w = 1.0 / (d2[i, ix] + 1e-6) * cw[ix]
        w = w / w.sum()
        frac[i] = (fr[ix] * w[:, None]).sum(axis=0)
        ok = np.isfinite(ksat[ix])
        if ok.any():
            lk, wk = np.log10(ksat[ix][ok]), w[ok] / w[ok].sum()
            o = np.argsort(lk)
            ks[i] = 10.0 ** np.interp(0.5, np.cumsum(wk[o]) - 0.5 * wk[o],
                                      lk[o])
    return frac, ks


def fit_predict(train, test, cols, use_cov):
    """One arm, one fold: class, fractions (both ways) and Ks.

    The classifier and the regressors take the covariates as extra features
    and handle a missing one themselves. The neighbour search takes only
    bulk density, and only for targets that have it, which is what the tool
    does (GshpReference(use_bd=True))."""
    c = list(cols) + (COVS if use_cov else [])
    xtr, xte = train[c].to_numpy(float), test[c].to_numpy(float)
    ok = np.isfinite(train[list(cols)].to_numpy(float)).all(axis=1)
    tr = train[ok]
    xtr = xtr[ok]
    clf = gbm().fit(xtr, tr["texture_class"].to_numpy())
    cls = clf.predict(xte)
    sand = gbm(regress=True).fit(xtr, tr["sand"].to_numpy()).predict(xte)
    clay = gbm(regress=True).fit(xtr, tr["clay"].to_numpy()).predict(xte)
    sand = np.clip(sand, 0, 100)
    clay = np.clip(clay, 0, 100 - sand)
    frac2 = np.column_stack([sand, 100.0 - sand - clay, clay])

    frac = np.full((len(test), 3), np.nan)
    ks = np.full(len(test), np.nan)
    has_bd = (np.isfinite(test["bd"].to_numpy(float)) if use_cov
              else np.zeros(len(test), bool))
    for kcols, sel in ((list(cols) + ["bd"], has_bd),
                       (list(cols), ~has_bd)):
        if not sel.any():
            continue
        ktr = tr[np.isfinite(tr[kcols].to_numpy(float)).all(axis=1)]
        a = ktr[kcols].to_numpy(float)
        b = test[kcols].to_numpy(float)[sel]
        mu, sd = a.mean(axis=0), a.std(axis=0) + 1e-12
        frac[sel], ks[sel] = knn((a - mu) / sd, (b - mu) / sd, ktr)
    return cls, frac, frac2, ks


def folds(tg, design):
    rng = np.random.default_rng(0)
    if design == "profile":
        p = np.array(sorted(tg.profile_id.unique()))
        rng.shuffle(p)
        return [("profile_id", set(x)) for x in np.array_split(p, N_FOLDS)]
    return [("source_db", {s}) for s in sorted(tg.source_db.unique())]


def score(name, truth, pred, base=None):
    ex = pred == truth
    gp = np.array([GROUP[a] == GROUP[b] for a, b in zip(pred, truth)])
    s = f"{name:8s} exact {ex.mean() * 100:5.1f}  group {gp.mean() * 100:5.1f}"
    if base is not None:
        b_ex, b_gp = base
        s += (f"   vs vg {(ex.mean() - b_ex.mean()) * 100:+5.1f} pp "
              f"(p={mcnemar(b_ex, ex):.3f})  group "
              f"{(gp.mean() - b_gp.mean()) * 100:+5.1f} (p={mcnemar(b_gp, gp):.3f})")
    print(s)
    return ex, gp


def ceiling(g, n_per_class, cap=300):
    """Cover & Hart bracket per predictor set, profile- and source-blocked.

    1NN in each (standardised) predictor space, on a reference capped per
    class -- the uniform prior the tool assumes -- so the brackets are
    comparable with verify_ceiling.py's.
    """
    from sklearn.model_selection import GroupKFold
    ref = pd.concat([x.sample(min(len(x), cap), random_state=0)
                     for _, x in g.groupby("texture_class")]).reset_index(drop=True)
    t = pd.concat([x.sample(min(len(x), n_per_class), random_state=1)
                   for _, x in ref.groupby("texture_class")])
    tpos = np.array([ref.index.get_loc(i) for i in t.index])
    truth = t.texture_class.to_numpy()
    print(f"reference {len(ref)} layers (<= {cap}/class), targets {len(t)}\n")
    for blocked in ("profile_id", "source_db"):
        print(f"--- 1NN blocked by {blocked} ---")
        print(f"{'predictors':24s} {'1NN':>6s}  {'Bayes accuracy':>18s}   "
              f"{'group 1NN':>9s}  {'Bayes group':>18s}")
        for arm, cols in ARMS.items():
            for use_cov in (False, True):
                c = list(cols) + (COVS if use_cov else [])
                x = ref[c].to_numpy(float)
                keep = np.isfinite(x).all(axis=1)
                mu, sd = x[keep].mean(0), x[keep].std(0) + 1e-12
                z = (x - mu) / sd
                pred = np.empty(len(t), object)
                for tr, te in GroupKFold(n_splits=N_FOLDS).split(
                        ref, groups=ref[blocked].to_numpy()):
                    sel = np.isin(tpos, te)
                    tr = tr[keep[tr]]
                    if not sel.any():
                        continue
                    a, b = z[tr], z[tpos[sel]]
                    d = ((b[:, None, :] - a[None, :, :]) ** 2).sum(2)
                    pred[sel] = ref.texture_class.to_numpy()[tr][d.argmin(1)]
                ok = np.array([p == q for p, q in zip(pred, truth)])
                gp = np.array([GROUP[p] == GROUP[q] for p, q in zip(pred, truth)])
                lo, hi = bayes_bracket(1 - ok.mean(), 12)
                glo, ghi = bayes_bracket(1 - gp.mean(), 4)
                name = arm + (" + depth, bd" if use_cov else "")
                print(f"{name:24s} {ok.mean() * 100:5.1f}%  "
                      f"[{lo * 100:5.1f}%, {hi * 100:5.1f}%]   "
                      f"{gp.mean() * 100:8.1f}%  [{glo * 100:5.1f}%, {ghi * 100:5.1f}%]")
        print()


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    measured_only = "--measured-only" in sys.argv
    design = args[0] if args else "profile"
    n_per_class = int(args[1]) if len(args) > 1 else 150
    if not os.path.exists(CACHE):
        raise SystemExit(f"{CACHE} not found -- run python ivg_fits.py first")

    g = st.load_reference_df("merged").reset_index(drop=True)
    g["source_db"] = g.source_db.fillna("KSSL")
    g["profile_id"] = g.profile_id.fillna("solo_" + g.layer_id.astype(str))
    cache = pd.read_csv(CACHE).set_index("layer_id")
    g = g.join(cache, on="layer_id")
    g = g[g.ivg_rmse.notna()]
    if measured_only:
        g = g[g.points == "measured"]
    g = g.reset_index(drop=True)
    counts = g.texture_class.value_counts()
    g["class_weight"] = (1.0 / counts)[g.texture_class].to_numpy()
    feats = features(g)
    g = pd.concat([g, feats], axis=1)

    if design == "ceiling":
        ceiling(g, n_per_class)
        return
    if design == "az":
        tg = g[g.source_db == "Babaeian_Arizona"].reset_index(drop=True)
        fold_list = [("source_db", {"Babaeian_Arizona"})]
    else:
        tg = pd.concat([x.sample(min(len(x), n_per_class), random_state=0)
                        for _, x in g.groupby("texture_class")]
                       ).reset_index(drop=True)
        fold_list = folds(tg, design)
    print(f"reference {len(g)} layers ({(g.points == 'measured').sum()} with "
          f"measured points); targets {len(tg)}; design {design}"
          f"{'; measured points only' if measured_only else ''}\n")

    n = len(tg)
    res = {(a, c): dict(cls=np.empty(n, object), frac=np.full((n, 3), np.nan),
                        frac2=np.full((n, 3), np.nan), ks=np.full(n, np.nan))
           for a in ARMS for c in (False, True)}
    for key, held in fold_list:
        sel = np.where(tg[key].isin(held))[0]
        if not len(sel):
            continue
        train = g[~g[key].isin(held)]
        test = tg.iloc[sel]
        for arm, cols in ARMS.items():
            for use_cov in (False, True):
                cls, frac, frac2, ks = fit_predict(train, test, cols, use_cov)
                r = res[(arm, use_cov)]
                r["cls"][sel] = cls
                r["frac"][sel] = frac
                r["frac2"][sel] = frac2
                r["ks"][sel] = ks
        print(f"  fold {sorted(held)[0]}: {len(sel)} targets", flush=True)

    truth = tg.texture_class.to_numpy()
    tfrac = tg[["sand", "silt", "clay"]].to_numpy()
    measured = (tg.points == "measured").to_numpy()
    dry = (tg.h_max_kpa > 1600).to_numpy()
    for use_cov in (False, True):
        print(f"\n=== texture class, {'curve + depth and bulk density' if use_cov else 'curve alone'} ===")
        base = None
        for arm in ARMS:
            base_now = score(arm, truth, res[(arm, use_cov)]["cls"], base)
            if arm == "vg":
                base = base_now
        print(f"--- targets with measured points (n={measured.sum()}) ---")
        b = None
        for arm in ARMS:
            x = score(arm, truth[measured], res[(arm, use_cov)]["cls"][measured], b)
            if arm == "vg":
                b = x
        if dry.sum() > 30:
            print(f"--- of those, dry end past 1500 kPa (n={dry.sum()}) ---")
            b = None
            for arm in ARMS:
                x = score(arm, truth[dry], res[(arm, use_cov)]["cls"][dry], b)
                if arm == "vg":
                    b = x

    print("\n=== particle fractions, MAE % (neighbour mean | sand+clay "
          "regressors, silt by difference) ===")
    print(f"{'arm':8s} {'cov':4s} " + "  ".join(f"{c:>21s}" for c in
                                                ("sand", "silt", "clay")))
    for arm in ARMS:
        for use_cov in (False, True):
            r = res[(arm, use_cov)]
            cells = []
            for j in range(3):
                a = np.nanmean(np.abs(r["frac"][:, j] - tfrac[:, j]))
                b = np.nanmean(np.abs(r["frac2"][:, j] - tfrac[:, j]))
                cells.append(f"{a:9.2f} | {b:9.2f}")
            print(f"{arm:8s} {'yes' if use_cov else 'no':4s} " +
                  "  ".join(cells))

    print("\n=== Ks, log10 RMSE (dex) and share within a factor of 10 ===")
    kt = tg.ksat_cmh.to_numpy()
    have = np.isfinite(kt)
    for arm in ARMS:
        line = f"{arm:8s}"
        for use_cov in (False, True):
            p = res[(arm, use_cov)]["ks"]
            ok = have & np.isfinite(p)
            e = np.log10(p[ok]) - np.log10(kt[ok])
            line += (f"   {'cov' if use_cov else 'curve':>5s} "
                     f"{np.sqrt(np.mean(e ** 2)):.3f} dex, "
                     f"{np.mean(np.abs(e) < 1) * 100:4.1f} % (n={ok.sum()})")
        print(line)


if __name__ == "__main__":
    main()
