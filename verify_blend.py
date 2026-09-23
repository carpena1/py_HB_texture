"""Can the classifier and the neighbour vote be combined by class?

The two halves of the hybrid are complementary rather than ranked. On
leave-one-soil-out (figures/validation_plots.py, v3) the classifier is better
on the common classes while the neighbour vote is much better on the rare
ones -- silt 59 % recall against 22 %, sandy clay 45 % against 22 % -- and,
unlike the usual rare-class trade, not at the cost of precision: silt F1 53
against 35. Overall the two are within 1.3 points.

Step 1 (this script, slow) caches BOTH probability vectors for every target
in each hold-out design, so that step 2 can compare combination rules
offline in seconds instead of re-running the pipeline for each one.

Step 2 (--analyse) compares, with the parameters always chosen on inner
folds and scored on an outer fold they never saw:

  gbm / knn     the two on their own, as they are today
  temperature   the classifier's probabilities scaled by one parameter
                first: it is overconfident, and blending a badly calibrated
                expert with a better one mostly measures the miscalibration
  linear        p = w p_gbm + (1-w) p_knn
  log pool      p proportional to p_gbm^a p_knn^(1-a)
  by frequency  a linear blend whose weight follows how well the reference
                covers the class, w_c = sigmoid(a + b log n_c): the
                hypothesis is "trust the neighbours where the reference is
                thin", stated as two parameters rather than a hand-picked
                list of rare classes
  by data       the same hypothesis inside the log pool, and with the right
                asymptote: the exponent on the classifier is
                a_c = 0.5 n_c / (n_c + k), so a class the reference barely
                covers leans on the neighbours and a well covered one tends
                to the plain geometric mean. One parameter, and k = 0 is
                exactly the geometric mean, so the fit can say whether the
                shrinkage earns anything

Judged on macro-F1 first -- a user reading a rare class cares about that,
not about overall accuracy -- then exact class and texture group.

With --heads both members use the fixed-head predictors
(swcc_texture.HEADS_CM) instead of the vG parameters, which says whether
that change and the ensemble add up.

Usage:  python verify_blend.py [layer|profile|source] [n_per_class] [n_mc]
                              [--heads]
        python verify_blend.py --analyse [design ...] [--heads]
"""

import os
import sys

import numpy as np
import pandas as pd

import swcc_texture as st
from verify_common import GROUP, mcnemar

H = np.arange(0.0, 150.0 + 0.1, 5.0)
N_FOLDS = 5
CLASSES = st.USDA_CLASSES
CACHE = os.path.join("figures", "cache", "blend_probs_{}{}.csv")


def vec(probs):
    """A class-probability dict as a vector over CLASSES."""
    return np.array([probs.get(c, 0.0) for c in CLASSES])


def compute(design, n_per_class, n_mc, heads=False):
    g = st.load_reference_df("merged").reset_index(drop=True)
    g["source_db"] = g.source_db.fillna("KSSL")
    g["profile_id"] = g.profile_id.fillna("solo_" + g.layer_id.astype(str))
    tg = pd.concat([x.sample(min(len(x), n_per_class), random_state=0)
                    for _, x in g.groupby("texture_class")]
                   ).reset_index(drop=True)
    rng = np.random.default_rng(0)
    if design == "layer":
        lay = np.array(tg.layer_id)
        rng.shuffle(lay)
        folds = [("layer_id", set(x)) for x in np.array_split(lay, N_FOLDS)]
    elif design == "profile":
        p = np.array(sorted(tg.profile_id.unique()))
        rng.shuffle(p)
        folds = [("profile_id", set(x)) for x in np.array_split(p, N_FOLDS)]
    else:
        folds = [("source_db", {s}) for s in sorted(tg.source_db.unique())]
    print(f"reference {len(g)} layers; targets {len(tg)}; design {design}; "
          f"n_mc={n_mc}", flush=True)

    pg = np.full((len(tg), len(CLASSES)), np.nan)
    pk = np.full((len(tg), len(CLASSES)), np.nan)
    for i, (key, held) in enumerate(folds, 1):
        sel = np.where(tg[key].isin(held))[0]
        if not len(sel):
            continue
        train = g[~g[key].isin(held)].reset_index(drop=True)
        fm = "heads" if heads else "vg"
        ref = st.GshpReference(df=train, use_bd=True, feature_mode=fm)
        clf = st.TextureGBM(df=train, covariates=["depth_cm", "bd"],
                            feature_mode=fm)
        for j in sel:
            r = tg.iloc[j]
            has_bd = bool(np.isfinite(r.bd))
            res = st.estimate(H, st.vg_theta(H, r.thetar, r.thetas,
                                             r.alpha_kpa, r.n),
                              ref=ref if has_bd else
                              st.GshpReference(df=train, feature_mode=fm),
                              clf=clf, n_mc=n_mc,
                              sample_type=r.sample_type, depth=r.depth_cm,
                              bulk_density=r.bd if has_bd else None)
            pg[j] = vec(res["class_probabilities"])
            pk[j] = vec(res["knn_class_probabilities"])
        print(f"  fold {i}/{len(folds)}: {len(sel)} targets", flush=True)

    out = pd.DataFrame({"layer_id": tg.layer_id, "source_db": tg.source_db,
                        "profile_id": tg.profile_id,
                        "texture_class": tg.texture_class})
    for k, c in enumerate(CLASSES):
        out[f"gbm_{c}"] = pg[:, k]
        out[f"knn_{c}"] = pk[:, k]
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    path = CACHE.format(design, "_heads" if heads else "")
    out.to_csv(path, index=False)
    print(f"wrote {path} ({len(out)} targets)")


# --------------------------------------------------------------------------
# step 2: combination rules, all fitted on inner folds only
# --------------------------------------------------------------------------

def norm(p):
    s = p.sum(axis=1, keepdims=True)
    return np.divide(p, s, out=np.full_like(p, 1.0 / p.shape[1]), where=s > 0)


def temper(p, t):
    """Temperature scaling: t > 1 softens an overconfident model."""
    with np.errstate(divide="ignore"):
        lp = np.log(np.clip(p, 1e-12, None)) / t
    e = np.exp(lp - lp.max(axis=1, keepdims=True))
    return norm(e)


def rules():
    """Each rule: name, parameter grid, and p(pg, pk, params) -> probabilities."""
    ws = np.round(np.arange(0.0, 1.01, 0.05), 2)
    ts = [0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
    yield "gbm", [()], lambda pg, pk, _: pg
    yield "knn", [()], lambda pg, pk, _: pk
    yield "temperature", [(t,) for t in ts], \
        lambda pg, pk, q: temper(pg, q[0])
    yield "linear", [(w,) for w in ws], \
        lambda pg, pk, q: norm(q[0] * pg + (1 - q[0]) * pk)
    yield "linear+temp", [(w, t) for w in ws for t in ts], \
        lambda pg, pk, q: norm(q[0] * temper(pg, q[1]) + (1 - q[0]) * pk)
    yield "log pool", [(a,) for a in ws], \
        lambda pg, pk, q: norm(np.clip(pg, 1e-12, None) ** q[0]
                               * np.clip(pk, 1e-12, None) ** (1 - q[0]))
    yield "by frequency", [(a, b) for a in (-4, -2, 0, 2, 4)
                           for b in (-2, -1, -0.5, 0, 0.5, 1, 2)], \
        lambda pg, pk, q: norm(_freq_w(q) * pg + (1 - _freq_w(q)) * pk)
    yield "by data", [(k,) for k in (0, 25, 50, 100, 200, 500, 1000, 2000,
                                     5000)], \
        lambda pg, pk, q: norm(np.clip(pg, 1e-12, None) ** _shrink(q[0])
                               * np.clip(pk, 1e-12, None) ** (1 - _shrink(q[0])))


_LOG_N = None
_COUNTS = None


def _shrink(k):
    """Exponent on the classifier, per class: 0.5 n/(n+k)."""
    return (0.5 * _COUNTS / (_COUNTS + k))[None, :]


def _freq_w(q):
    a, b = q
    return 1.0 / (1.0 + np.exp(-(a + b * _LOG_N)))[None, :]


def macro_f1(truth, pred):
    f = []
    for c in CLASSES:
        tp = np.sum((pred == c) & (truth == c))
        d = np.sum(pred == c) + np.sum(truth == c)
        f.append(2 * tp / d if d else 0.0)
    return float(np.mean(f)) * 100


def score(truth, pred):
    ex = (pred == truth)
    gp = np.array([GROUP[a] == GROUP[b] for a, b in zip(pred, truth)])
    return ex, gp, macro_f1(truth, pred)


def analyse(design, heads=False):
    global _LOG_N
    path = CACHE.format(design, "_heads" if heads else "")
    if not os.path.exists(path):
        print(f"  {path} not built yet", file=sys.stderr)
        return
    d = pd.read_csv(path)
    pg = d[[f"gbm_{c}" for c in CLASSES]].to_numpy()
    pk = d[[f"knn_{c}" for c in CLASSES]].to_numpy()
    truth = d.texture_class.to_numpy()
    global _COUNTS
    counts = st.load_reference_df("merged").texture_class.value_counts()
    _COUNTS = np.array([counts.get(c, 1) for c in CLASSES], float)
    _LOG_N = np.log(np.array([counts.get(c, 1) for c in CLASSES], float))
    _LOG_N = (_LOG_N - _LOG_N.mean()) / _LOG_N.std()

    # Outer folds by the same unit the design hides, so a rule is never
    # tuned on the group it is scored on.
    key = {"layer": "layer_id", "profile": "profile_id",
           "source": "source_db"}[design]
    units = np.array(sorted(d[key].unique()))
    rng = np.random.default_rng(0)
    rng.shuffle(units)
    outer = [set(u) for u in np.array_split(units, N_FOLDS)]

    print(f"\n=== {design} design, {len(d)} targets, "
          f"{'fixed heads' if heads else 'vG parameters'} "
          f"(parameters fitted on inner folds only) ===")
    print(f"{'rule':14s} {'macro-F1':>9s} {'exact':>7s} {'group':>7s}   "
          f"{'vs hybrid today':>28s}   chosen")
    base_ex = base_f1 = None
    for name, grid, f in rules():
        pred = np.empty(len(d), object)
        chosen = []
        for held in outer:
            te = d[key].isin(held).to_numpy()
            tr = ~te
            if not te.any() or not tr.any():
                continue
            best, best_f1 = grid[0], -1.0
            for q in grid:                      # fit on the outer-fold complement
                pi = f(pg[tr], pk[tr], q)
                f1 = macro_f1(truth[tr],
                              np.array(CLASSES)[pi.argmax(axis=1)])
                if f1 > best_f1:
                    best, best_f1 = q, f1
            pred[te] = np.array(CLASSES)[f(pg[te], pk[te], best).argmax(axis=1)]
            chosen.append(best)
        ex, gp, f1 = score(truth, pred)
        if name == "gbm":
            base_ex, base_gp, base_f1 = ex, gp, f1
        cmp = ""
        if base_ex is not None and name != "gbm":
            cmp = (f"F1 {f1 - base_f1:+5.1f}  exact {(ex.mean() - base_ex.mean()) * 100:+5.1f} pp "
                   f"(p={mcnemar(base_ex, ex):.3f})")
        uniq = sorted({tuple(round(float(x), 2) for x in np.atleast_1d(c))
                       for c in chosen})
        shown = ", ".join(str(u[0] if len(u) == 1 else u) for u in uniq[:3])
        print(f"{name:14s} {f1:9.1f} {ex.mean() * 100:6.1f} % "
              f"{gp.mean() * 100:6.1f} %   {cmp:>28s}   "
              f"{shown}{'...' if len(uniq) > 3 else ''}")


def main():
    heads = "--heads" in sys.argv
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if "--analyse" in sys.argv:
        for design in (args or ["layer", "profile", "source"]):
            analyse(design, heads)
        return
    design = args[0] if args else "layer"
    compute(design, int(args[1]) if len(args) > 1 else 150,
            int(args[2]) if len(args) > 2 else 40, heads)


if __name__ == "__main__":
    main()
