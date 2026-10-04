"""What a new sensor network can expect, and what improves it.

Three questions behind the definition of a "source" in the hold-out tests.
Every prediction runs as a network would use the tool: depth and bulk density
supplied (TextureGBM covariates, use_bd reference), curves regenerated from
each target's vG parameters and refitted (st.estimate, n_mc = 20).

  a  proximity   With the target's laboratory hidden, are soils predicted
                 better when reference soils from OTHER laboratories lie
                 close by? Distance to the nearest such soil; texture
                 correctness and Ks error regressed on log distance with
                 texture-class fixed effects; bootstrap over sources.
  b  learning    A network adds lab-verified sites of its own: for each
                 source with enough profiles, fixed test profiles (one layer
                 each) are predicted with 0, 1, 2, 5, 10, 20, 50 or all of its
                 other profiles added back (three random draws per k).
  c  units       Compilations split into their studies where the metadata
                 allow (ETH literature by reference, UNSODA by method, the
                 Australian database by reference, KSSL by its two 33-kPa
                 method eras, EU-HYDI contributors by data provider, AfSPDB
                 by country and survey): hiding the target's own
                 study only, against hiding the whole compilation.

  d  weighting   The same learning curve with the network's own verified
                 profiles made to count: preferred in the neighbour search
                 (other networks' soils lambda = 1 or 3 standard units further
                 away), weighted 20 times in the classifier, both, or a local
                 Ks offset from the verified cores' residuals (shrunk by
                 n / (n + 5)) applied to the new-network prediction.

Usage:  python verify_network.py [a|b|c|d ...] [--reuse]   (default: all four)
Per-target predictions go to figures/cache/ (git-ignored: EU-HYDI layers).
Folds run in parallel (verify_common.map_folds).
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.neighbors import BallTree

import swcc_texture as st
from verify_common import GROUP, map_folds

H = np.arange(0.0, 150.0 + 0.1, 5.0)
N_MC = 20
LOG2 = np.log10(2.0)
CACHE = os.path.join("figures", "cache")
GSHP_RAW = os.path.join("data", "WRC_dataset_surya_et_al_2021_final.csv")
KSSL_DB = os.path.join("data", "kssl_raw", "NCSSLabDataMartSQLite.sqlite3")
K_GRID = [0, 1, 2, 5, 10, 20, 50]
N_TEST, MIN_POOL, REPS = 40, 20, 3


def load():
    g = st.load_reference_df(st.DEFAULT_REFERENCE).reset_index(drop=True)
    g["source_db"] = g.source_db.fillna("KSSL")
    g["profile_id"] = g.profile_id.fillna("solo_" + g.layer_id.astype(str))
    g.loc[~g.bd.between(0.1, 2.3), "bd"] = np.nan
    # study: the finest unit the metadata support
    raw = pd.read_csv(GSHP_RAW, low_memory=False, encoding="latin-1",
                      usecols=["layer_id", "reference", "method"]).drop_duplicates("layer_id")
    g = g.merge(raw, on="layer_id", how="left")
    g["study"] = g.source_db
    for src, col in (("ETH_Literature", "reference"), ("UNSODA", "method"),
                     ("Australian_database", "reference")):
        m = g.source_db == src
        g.loc[m, "study"] = src + ":" + g.loc[m, col].fillna("?").astype(str)
    # The NCSS Characterization Database is KSSL's with cooperating
    # laboratories, but its method codes are KSSL handbook procedures that the
    # cooperators share, so they do not name the laboratory. The retained
    # layers fall into two 33-kPa procedures (4B1c and DbWR1, matching the two
    # sample-number series): two method eras, kept as two units.
    # EU-HYDI contributors that are themselves compilations (the HYPRES
    # sets) are split by data provider (CONTACT_P), or country where the
    # provider is not given; AfSPDB by country and survey, from its profile
    # identifiers (e.g. "KE SOTER_...", "KE wasp_...", "NG ...").
    m = g.source_db.str.startswith("EUHYDI_")
    if m.any():
        try:
            import prepare_euhydi as pe
            gen = pe.table("GENERAL")[["PROFILE_ID", "CONTACT_P", "ISO_COUNTRY"]]
            prov = gen.CONTACT_P.where(~gen.CONTACT_P.astype(str).str.strip()
                                       .isin(["NAME", "ND", "nan", ""]), gen.ISO_COUNTRY)
            prov.index = "EUHYDI_" + gen.PROFILE_ID.astype(str)
            prov = prov[~prov.index.duplicated()]
            g.loc[m, "study"] = (g.loc[m, "source_db"] + ":"
                                 + g.loc[m, "profile_id"].map(prov).fillna("?").astype(str))
        except Exception as e:                       # raw database not here
            print(f"  EU-HYDI providers not available ({e}); contributors kept whole",
                  file=sys.stderr)
    m = g.source_db == "AfSPDB"
    tok = g.loc[m, "profile_id"].astype(str).str.extract(r"^([A-Z]{2})\s*([A-Za-z]*)")
    g.loc[m, "study"] = "AfSPDB:" + tok[0].fillna("?") + " " + tok[1].fillna("")
    if os.path.exists(KSSL_DB):
        import sqlite3
        with sqlite3.connect(KSSL_DB) as con:
            mth = pd.read_sql("SELECT layer_key, water_retention_thirdbar_method AS m "
                              "FROM lab_physical_properties_vw", con)
        mth = (mth.dropna().assign(layer_id=lambda x: "KSSL_" + x.layer_key.astype(str))
               .drop_duplicates("layer_id").set_index("layer_id").m.astype(str).str.strip())
        m = g.source_db == "KSSL"
        g.loc[m, "study"] = "KSSL:" + g.loc[m, "layer_id"].map(mth).fillna("?")
    return g


def predict(train, tg, idx):
    train = train.reset_index(drop=True)
    ref = st.GshpReference(df=train, use_bd=True)
    clf = st.TextureGBM(df=train, covariates=["depth_cm", "bd"])
    out = []
    for j in idx:
        r = tg.loc[j]
        res = st.estimate(H, st.vg_theta(H, r.thetar, r.thetas, r.alpha_kpa, r.n),
                          ref=ref, n_mc=N_MC, clf=clf, sample_type=r.sample_type,
                          depth=r.depth_cm, bulk_density=r.bd)
        out.append(dict(j=j, cls=res["texture_class"], ks=res["ksat"]["median_cmh"]))
    return out


def hide_fold(item, g, tg):
    unit, value = item
    idx = tg.index[tg[unit] == value]
    out = predict(g[g[unit] != value], tg, idx)
    for o in out:
        o["unit"] = unit
    return out


def correct(d):
    return (d.cls == d.texture_class).to_numpy(), np.array(
        [GROUP[a] == GROUP[b] for a, b in zip(d.cls, d.texture_class)])


def ks_err(d):
    k = d[d.ksat_cmh > 0]
    return np.abs(np.log10(k.ks.to_numpy(float)) - np.log10(k.ksat_cmh.to_numpy(float))), k


# ------------------------------------------------------------------ a / c
def targets(g):
    pool = g[g.depth_cm.notna() & g.bd.notna()]
    return pd.concat([x.sample(min(len(x), 150), random_state=0)
                      for _, x in pool.groupby("texture_class")])


def hidden_runs(g, tg):
    path = os.path.join(CACHE, "network_hidden.csv")
    if os.path.exists(path) and "--reuse" in sys.argv:
        return pd.read_csv(path, index_col=0)
    items = [("source_db", s) for s in sorted(tg.source_db.unique())]
    split = tg.study != tg.source_db
    items += [("study", s) for s in sorted(tg.study[split].unique())]
    est = pd.DataFrame([r for part in map_folds(hide_fold, items, g, tg) for r in part])
    est.to_csv(path)
    return est


def boot_coef(fit, d, n=1000, seed=1):
    rng = np.random.default_rng(seed)
    srcs = d.source_db.unique()
    by = {s: np.where(d.source_db.to_numpy() == s)[0] for s in srcs}
    vals = []
    for _ in range(n):
        ii = np.concatenate([by[s] for s in rng.choice(srcs, len(srcs))])
        try:
            vals.append(fit(d.iloc[ii]))
        except ValueError:
            pass
    return np.percentile(vals, [2.5, 97.5]), np.mean(np.array(vals) <= 0)


DIST_BINS = [0, 10, 50, 200, 1000, 1e5]
DIST_LABELS = ["< 10 km", "10-50 km", "50-200 km", "200-1000 km", "> 1000 km"]


def proximity_frame(g, tg, est):
    """Laboratory-hidden predictions with the distance to the nearest
    reference soil from another laboratory: (all soils, soils with Ks)."""
    d = est[est.unit == "source_db"].set_index("j").join(tg)
    d = d[d.lat.notna()]
    dist = np.full(len(d), np.nan)
    geo = g[g.lat.notna()]
    for s, x in d.groupby("source_db"):
        other = geo[geo.source_db != s]
        tree = BallTree(np.radians(other[["lat", "lon"]].to_numpy()), metric="haversine")
        dd, _ = tree.query(np.radians(x[["lat", "lon"]].to_numpy()), k=1)
        dist[d.index.get_indexer(x.index)] = dd[:, 0] * 6371.0
    d = d.assign(km=np.maximum(dist, 1.0))
    ex, gp = correct(d)
    d = d.assign(ex=ex, gp=gp)
    d["bin"] = pd.cut(d.km, DIST_BINS, labels=DIST_LABELS)
    e, k = ks_err(d)
    return d, k.assign(e=e)


def test_a(g, tg, est):
    d, k = proximity_frame(g, tg, est)
    print(f"\n=== (a) PROXIMITY, laboratory hidden: {len(d)} soils with coordinates, "
          f"{d.source_db.nunique()} sources ===")
    print("  nearest other-laboratory soil    soils   exact   group   | Ks soils  within 10x  median factor")
    for b, x in d.groupby("bin", observed=True):
        y = k[k.bin == b]
        print(f"  {b:30s} {len(x):6d}  {x.ex.mean() * 100:5.1f} %  {x.gp.mean() * 100:5.1f} %  | "
              f"{len(y):6d}   {np.mean(y.e <= 1) * 100:5.1f} %      x{10 ** np.median(y.e):.1f}"
              if len(y) else f"  {b:30s} {len(x):6d}  {x.ex.mean() * 100:5.1f} %  {x.gp.mean() * 100:5.1f} %  |")

    def slope_cls(dd, y):
        X = np.column_stack([np.log10(dd.km), pd.get_dummies(dd.texture_class).to_numpy(float)])
        return LogisticRegression(C=1e6, max_iter=2000).fit(X, dd[y]).coef_[0][0]

    for y, lab in (("ex", "exact class"), ("gp", "texture group")):
        b = slope_cls(d, y)
        (lo, hi), p = boot_coef(lambda dd: slope_cls(dd, y), d)
        print(f"  {lab}: log-odds per tenfold distance {b:+.3f} [95 % CI {lo:+.3f}, {hi:+.3f}]"
              f" (negative = closer is better; share of resamples >= 0: {1 - p:.3f})")

    def slope_ks(dd):
        X = np.column_stack([np.log10(dd.km), pd.get_dummies(dd.texture_class).to_numpy(float)])
        return LinearRegression().fit(X, dd.e).coef_[0]
    b = slope_ks(k)
    (lo, hi), p = boot_coef(slope_ks, k)
    print(f"  Ks |log10 error| per tenfold distance {b:+.3f} dex [95 % CI {lo:+.3f}, {hi:+.3f}]"
          f" (positive = closer is better)")
    print(f"  rank correlation, distance vs Ks error: {spearmanr(k.km, k.e).statistic:+.3f}")


def test_c(tg, est):
    print("\n=== (c) UNITS: own study hidden vs whole compilation hidden ===")
    a = est[est.unit == "source_db"].set_index("j")
    b = est[est.unit == "study"].set_index("j")
    common = a.index.intersection(b.index)
    d = tg.loc[common].assign(cls_a=a.loc[common].cls, cls_b=b.loc[common].cls,
                              ks_a=a.loc[common].ks, ks_b=b.loc[common].ks)
    from verify_common import mcnemar
    for src, x in [("all split compilations", d)] + list(d.groupby("source_db")):
        ea, eb = (x.cls_a == x.texture_class).to_numpy(), (x.cls_b == x.texture_class).to_numpy()
        ga = np.array([GROUP[p] == GROUP[t] for p, t in zip(x.cls_a, x.texture_class)])
        gb = np.array([GROUP[p] == GROUP[t] for p, t in zip(x.cls_b, x.texture_class)])
        k = x[x.ksat_cmh > 0]
        line = (f"  {src:24s} n={len(x):4d}  exact {ea.mean() * 100:5.1f} -> {eb.mean() * 100:5.1f} %"
                f" (p={mcnemar(ea, eb):.2f})  group {ga.mean() * 100:5.1f} -> {gb.mean() * 100:5.1f} %")
        if len(k) > 5:
            e_a = np.abs(np.log10(k.ks_a / k.ksat_cmh)); e_b = np.abs(np.log10(k.ks_b / k.ksat_cmh))
            line += f"  Ks within 10x {np.mean(e_a <= 1) * 100:4.0f} -> {np.mean(e_b <= 1) * 100:4.0f} %"
        print(line)
    print(f"  ({tg.study[tg.study != tg.source_db].nunique()} studies or methods inside "
          f"{tg.source_db[tg.study != tg.source_db].nunique()} compilations; WOSIS, HYBRAS and "
          f"Russia_EGRPR carry no finer reference and stay whole)")


# ------------------------------------------------------------------ b
def added_profiles(pool, k, rep):
    """The k verified profiles a network adds (the same draw in every test)."""
    if k == "all":
        return set(pool)
    rng = np.random.default_rng(1000 * rep + k)
    return set(rng.choice(pool, k, replace=False)) if k else set()


def learn_fold(item, g, plan):
    src, k, rep = item
    test_layers, pool = plan[src]
    added = added_profiles(pool, k, rep)
    train = g[(g.source_db != src) | g.profile_id.isin(added)]
    tg = g.set_index("layer_id").loc[test_layers].reset_index()
    out = predict(train, tg, tg.index)
    for o, lid in zip(out, test_layers):
        o.update(source=src, k=k, rep=rep, layer_id=lid)
    return out


def make_plan(g):
    """Per source with enough profiles: N_TEST fixed test profiles (one layer
    each, carrying depth and bulk density) and the pool of its other profiles."""
    rng = np.random.default_rng(0)
    ok = g[g.depth_cm.notna() & g.bd.notna()]
    plan = {}
    for s, x in ok.groupby("source_db"):
        profs = np.array(sorted(x.profile_id.unique()))
        if len(profs) < N_TEST + MIN_POOL:
            continue
        rng.shuffle(profs)
        test_p = profs[:N_TEST]
        layers = [x[x.profile_id == p].sample(1, random_state=0).layer_id.iloc[0] for p in test_p]
        pool = list(g[(g.source_db == s) & ~g.profile_id.isin(test_p)].profile_id.unique())
        plan[s] = (layers, pool)
    return plan


def test_b(g):
    print("\n=== (b) LEARNING CURVE: the network's own verified profiles added ===")
    plan = make_plan(g)
    items = []
    for s, (_, pool) in plan.items():
        items.append((s, 0, 0))
        items += [(s, k, r) for k in K_GRID[1:] if k < len(pool) for r in range(REPS)]
        items.append((s, "all", 0))
    path = os.path.join(CACHE, "network_learning.csv")
    if os.path.exists(path) and "--reuse" in sys.argv:
        est = pd.read_csv(path)
    else:
        est = pd.DataFrame([r for part in map_folds(learn_fold, items, g, plan) for r in part])
        est.to_csv(path, index=False)
    est["k"] = est.k.astype(str)
    truth = g.set_index("layer_id")
    est = est.join(truth[["texture_class", "ksat_cmh"]], on="layer_id")
    print(f"  {len(plan)} sources, {N_TEST} test profiles each (one layer per profile); "
          f"k = verified profiles of the same source added")
    ex, gp = correct(est)
    est = est.assign(ex=ex, gp=gp,
                     e=np.where(est.ksat_cmh > 0, np.abs(np.log10(est.ks / est.ksat_cmh)), np.nan))
    # one value per source and k (mean over draws), then pooled with a bootstrap over sources
    per = est.groupby(["source", "k"]).agg(ex=("ex", "mean"), gp=("gp", "mean"),
                                           w10=("e", lambda v: np.mean(v.dropna() <= 1) if v.notna().sum() else np.nan),
                                           w2=("e", lambda v: np.mean(v.dropna() <= LOG2) if v.notna().sum() else np.nan),
                                           fac=("e", lambda v: 10 ** np.nanmedian(v) if v.notna().sum() else np.nan)).reset_index()
    base = per[per.k == "0"].set_index("source")
    order = [str(k) for k in K_GRID] + ["all"]
    rb = np.random.default_rng(2)
    print("  k     sources   exact (gain vs k=0, 95 % CI)        group (gain)          Ks within 10x  within 2x  factor")
    for k in order:
        x = per[per.k == k].set_index("source")
        if not len(x):
            continue
        s = x.index
        gain = lambda c: (x[c] - base.loc[s, c]).to_numpy()
        def ci(c):
            v = gain(c); v = v[np.isfinite(v)]
            bs = [np.mean(rb.choice(v, len(v))) for _ in range(2000)]
            return np.mean(v), np.percentile(bs, 2.5), np.percentile(bs, 97.5)
        ge, gg = ci("ex"), ci("gp")
        print(f"  {k:5s} {len(x):5d}    {x.ex.mean() * 100:5.1f} % ({ge[0] * 100:+.1f} [{ge[1] * 100:+.1f}, {ge[2] * 100:+.1f}])"
              f"   {x.gp.mean() * 100:5.1f} % ({gg[0] * 100:+.1f} [{gg[1] * 100:+.1f}, {gg[2] * 100:+.1f}])"
              f"   {np.nanmean(x.w10) * 100:5.1f} %       {np.nanmean(x.w2) * 100:5.1f} %   x{np.nanmedian(x.fac):.1f}")



# ------------------------------------------------------------------ d
LAMBDAS = (1.0, 3.0)      # extra standard units between a target and other networks' soils
DUP = 20                  # classifier weight of the network's own rows
SHRINK = 5                # Ks offset shrunk by n / (n + SHRINK) verified cores


def _neighbors_with_prior(self, thetar, thetas, alpha, n, k, depth=None, om=None,
                          m=None, sample_type=None, bd=None, andic=None,
                          require_ksat=False):
    """st.GshpReference.neighbors with one addition: self._prior (squared
    standard units per row, or None) is added to every distance, so the
    network's own rows can be preferred. Test-only; the tool is unchanged."""
    f = list(st._head_features(thetar, thetas, alpha, n, st.HEADS_KPA).ravel())
    if self.use_bd:
        f.append(np.log10(1.0 + bd))
    f = (np.array(f) - self.mean) / self.std
    d2 = ((self.z - f) ** 2).sum(axis=1)
    if andic is not None:
        d2 = d2 + (st.ANDIC_LAMBDA * (self.andic - andic)) ** 2
    prior = getattr(self, "_prior", None)
    if prior is not None:
        d2 = d2 + prior
    d = np.sqrt(d2)
    excluded = getattr(self, "_excluded", None)
    if excluded is not None:
        d = np.where(excluded, np.inf, d)
    if sample_type is not None:
        d = np.where(self.sample_type == sample_type, d, np.inf)
    if require_ksat:
        d = np.where(self.has_ksat, d, np.inf)
        k = min(k, int(np.isfinite(d).sum()))
        if k == 0:
            return np.array([], dtype=int), np.array([])
    idx = np.argpartition(d, k)[:k] if k < len(d) else np.argsort(d)[:k]
    w = 1.0 / (d[idx] ** 2 + 1e-6)
    return idx, w / w.sum()


def weight_fold(item, g, plan):
    """One (source, k, draw): the base tool and the weighted variants, on the
    same test layers. k = 0 is the new-network baseline."""
    st.GshpReference.neighbors = _neighbors_with_prior
    src, k, rep = item
    test_layers, pool = plan[src]
    added = added_profiles(pool, k, rep)
    train = g[(g.source_db != src) | g.profile_id.isin(added)].reset_index(drop=True)
    own = train.profile_id.isin(added).to_numpy() & (train.source_db == src).to_numpy()
    ref = st.GshpReference(df=train, use_bd=True)
    # the bulk-density reference drops rows without a density: align by layer
    own_ref = np.isin(ref.layer_id, train.layer_id[own].to_numpy())
    clf0 = st.TextureGBM(df=train, covariates=["depth_cm", "bd"])
    arms = [("base", None, clf0)]
    if k:
        clfw = st.TextureGBM(df=pd.concat([train] + [train[own]] * (DUP - 1), ignore_index=True),
                             covariates=["depth_cm", "bd"])
        arms += [(f"knn x{lam:g}", lam, clf0) for lam in LAMBDAS]
        arms += [(f"gbm x{DUP}", None, clfw), (f"both", LAMBDAS[-1], clfw)]
    tg = g.set_index("layer_id").loc[test_layers].reset_index()
    out = []
    for name, lam, clf in arms:
        ref._prior = None if lam is None else np.where(own_ref, 0.0, lam ** 2)
        for j in tg.index:
            r = tg.loc[j]
            res = st.estimate(H, st.vg_theta(H, r.thetar, r.thetas, r.alpha_kpa, r.n),
                              ref=ref, n_mc=N_MC, clf=clf, sample_type=r.sample_type,
                              depth=r.depth_cm, bulk_density=r.bd)
            out.append(dict(source=src, k=k, rep=rep, arm=name, layer_id=r.layer_id,
                            cls=res["texture_class"],
                            knn=next(iter(res["knn_class_probabilities"])),
                            ks=res["ksat"]["median_cmh"]))
    return out


def pool_fold(src, g, plan):
    """Ks of one core per pool profile, predicted as for a new network: the
    residuals the verified cores would give for a local Ks offset."""
    _, pool = plan[src]
    p = g[(g.source_db == src) & g.profile_id.isin(pool) & (g.ksat_cmh > 0)]
    p = p.groupby("profile_id").head(1)
    if not len(p):
        return []
    tg = p.reset_index(drop=True)
    out = predict(g[g.source_db != src], tg, tg.index)
    return [dict(source=src, profile_id=tg.profile_id[o["j"]],
                 resid=np.log10(tg.ksat_cmh[o["j"]] / o["ks"])) for o in out]


def test_d(g):
    print("\n=== (d) OWN-SITE WEIGHTING: a network's few verified profiles made to count ===")
    plan = make_plan(g)
    path = os.path.join(CACHE, "network_weighting.csv")
    ppath = os.path.join(CACHE, "network_pool_resid.csv")
    if os.path.exists(path) and "--reuse" in sys.argv:
        est, res = pd.read_csv(path), pd.read_csv(ppath)
    else:
        items = []
        for s, (_, pool) in plan.items():
            items.append((s, 0, 0))
            items += [(s, k, r) for k in K_GRID[1:] if k < len(pool) for r in range(REPS)]
        est = pd.DataFrame([r for part in map_folds(weight_fold, items, g, plan) for r in part])
        res = pd.DataFrame([r for part in map_folds(pool_fold, list(plan), g, plan) for r in part])
        est.to_csv(path, index=False)
        res.to_csv(ppath, index=False)
    per = weighting_per_source(g, plan, est, res)
    b0 = per[(per.k == 0)].set_index("source")
    rb = np.random.default_rng(3)

    def gain(x, c):
        v = (x[c] - b0.loc[x.index, c]).to_numpy()
        v = v[np.isfinite(v)]
        if not len(v):
            return "      —       "
        bs = [np.mean(rb.choice(v, len(v))) for _ in range(2000)]
        lo, hi = np.percentile(bs, [2.5, 97.5])
        sig = "*" if lo > 0 or hi < 0 else " "
        sc = 100 if c != "lerr" else 1
        fmt = "+.1f" if c != "lerr" else "+.2f"
        return f"{np.mean(v) * sc:{fmt}}{sig}[{lo * sc:{fmt}},{hi * sc:{fmt}}]"
    print(f"  {len(plan)} sources; gains are against the new network (k = 0), paired, 95 % CI over sources;"
          f" * = interval excludes 0")
    k0 = per[per.k == 0]
    print(f"  k = 0: exact {k0.ex.mean() * 100:.1f} %, group {k0.gp.mean() * 100:.1f} %, Ks within 10x "
          f"{np.nanmean(k0.w10) * 100:.1f} %, within 2x {np.nanmean(k0.w2) * 100:.1f} %")
    arms = ["base", "knn x1", "knn x3", f"gbm x{DUP}", "both", "Ks offset"]
    print(f"  {'k':>3s} {'arm':10s} {'exact gain':22s} {'group gain':22s} {'Ks within 10x':22s} {'Ks within 2x':22s} {'Ks log err':18s}")
    for k in K_GRID[1:]:
        for a in arms:
            x = per[(per.k == k) & (per.arm == a)].set_index("source")
            if not len(x):
                continue
            print(f"  {k:3d} {a:10s} {gain(x, 'ex'):22s} {gain(x, 'gp'):22s} {gain(x, 'w10'):22s} "
                  f"{gain(x, 'w2'):22s} {gain(x, 'lerr'):18s}")


def weighting_per_source(g, plan, est, res):
    """Per source, k and arm: texture and Ks scores, the local Ks offset arm
    built from the verified cores' residuals (one value per source: the mean
    over the draws)."""
    truth = g.set_index("layer_id")[["texture_class", "ksat_cmh"]]
    est = est.join(truth, on="layer_id")
    # local Ks offset: k = 0 predictions times 10^offset from the k verified cores
    base0 = est[(est.k == 0) & (est.arm == "base")]
    rows = []
    for (s, k, rep), _ in est[est.k > 0].groupby(["source", "k", "rep"]):
        added = added_profiles(plan[s][1], k, rep)
        r = res[(res.source == s) & res.profile_id.isin(added)].resid
        off = np.median(r) * len(r) / (len(r) + SHRINK) if len(r) else 0.0
        x = base0[base0.source == s].copy()
        x["k"], x["rep"], x["arm"], x["ks"] = k, rep, "Ks offset", x.ks * 10 ** off
        rows.append(x)
    est = pd.concat([est] + rows, ignore_index=True)
    est["ex"] = est.cls == est.texture_class
    est["gp"] = [GROUP[a] == GROUP[b] for a, b in zip(est.cls, est.texture_class)]
    est["kx"] = est.knn == est.texture_class
    est["e"] = np.where(est.ksat_cmh > 0, np.abs(np.log10(est.ks / est.ksat_cmh)), np.nan)
    per = est.groupby(["source", "k", "arm"]).agg(
        ex=("ex", "mean"), gp=("gp", "mean"), kx=("kx", "mean"),
        w10=("e", lambda v: np.mean(v.dropna() <= 1) if v.notna().sum() else np.nan),
        w2=("e", lambda v: np.mean(v.dropna() <= LOG2) if v.notna().sum() else np.nan),
        lerr=("e", lambda v: np.nanmedian(v) if v.notna().sum() else np.nan)).reset_index()
    return per


def main():
    parts = [a for a in sys.argv[1:] if a in ("a", "b", "c", "d")] or ["a", "b", "c", "d"]
    os.makedirs(CACHE, exist_ok=True)
    g = load()
    if {"a", "c"} & set(parts):
        tg = targets(g)
        est = hidden_runs(g, tg)
        if "a" in parts:
            test_a(g, tg, est)
        if "c" in parts:
            test_c(tg, est)
    if "b" in parts:
        test_b(g)
    if "d" in parts:
        test_d(g)


if __name__ == "__main__":
    main()
