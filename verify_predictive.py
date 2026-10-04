"""Is the tool predictive for a new sensor network, and how accurate is it?

The use case: an in-situ sensor network delivers retention curves from a
protocol the reference has never seen; depth and usually bulk density are
known from sampling at installation, and the cores taken then can be
analysed for texture and Ks to serve as ground truth. This script answers, with statistics:

  * what accuracy to expect for texture (class, class in the top two, group,
    fractions) and Ks (within 2x and 10x, typical factor, ranking, band
    coverage), for the neighbour estimate and the physical one;
  * whether the curve is predictive at all: every measure is compared, paired
    on the same soils, with what a user already knows without it -- depth and
    bulk density alone (a gradient-boosted model on those two), chance, and the
    reference's median Ks;
  * how far that accuracy is from the best any method could do (texture: the
    Cover & Hart ceiling, verify_ceiling.py; Ks: the spread between replicate
    cores of the same soil in the reference);
  * how much it varies from one source to the next, which is what a single
    new network will see;
  * whether the tool's own confidence signals (class probability, the
    neighbour vote agreeing, the physical Ks agreeing) pick out the
    predictions to trust.

A "source" here is one study or laboratory -- one measurement protocol --
not a database: compilations (ETH literature, UNSODA, the Australian set)
are split into their studies or methods where the metadata allow
(verify_network.load). A new in-situ network brings a protocol the
reference does not contain (every reference curve is a laboratory curve),
so hiding the target's study is the honest test for it, and possibly an
optimistic one until in-situ curves with cores are available.

Scenarios (targets: 150 per class from the soils that carry depth and bulk
density, curves regenerated from their vG parameters and refitted):
  source   the target's study hidden: a new network
  profile  only the target's profile hidden: an established network, whose
           own lab-verified sites (cores analysed for texture and Ks) are in
           the reference in quantity (add_local_data.py). verify_network.py
           (b, d) measures the stages between: a handful of verified sites
           changes nothing; about 20 begin to help and 50 clearly do.
Uncertainty: 95 % intervals from a bootstrap over studies (2,000
resamples), the study being the independent unit; differences from the
baselines are bootstrapped the same way, paired.

Usage:  python verify_predictive.py [n_per_class] [n_mc] [--reuse]
Per-target predictions go to figures/cache/ (git-ignored: EU-HYDI layers).
"""

import itertools
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor

import swcc_texture as st
from verify_common import GROUP, map_folds

H = np.arange(0.0, 150.0 + 0.1, 5.0)
CACHE = os.path.join("figures", "cache", "predictive.csv")
N_BOOT = 2000
LOG2 = np.log10(2.0)


def fold(item, g, tg, n_mc):
    design, col, held = item
    train = g[~g[col].isin(held)].reset_index(drop=True)
    ref_bd = st.GshpReference(df=train, use_bd=True)
    clf_cov = st.TextureGBM(df=train, covariates=["depth_cm", "bd"])
    # What a user knows without the curve: depth and bulk density.
    cov = train[train.depth_cm.notna() & train.bd.notna()]
    xc = np.column_stack([np.log10(1 + cov.depth_cm), cov.bd])
    base_cls = HistGradientBoostingClassifier(
        max_iter=300, learning_rate=0.06, early_stopping=True, random_state=0,
        class_weight="balanced").fit(xc, cov.texture_class)
    ck = cov[cov.ksat_cmh > 0]
    base_ks = HistGradientBoostingRegressor(
        max_iter=300, learning_rate=0.06, early_stopping=True, random_state=0
    ).fit(np.column_stack([np.log10(1 + ck.depth_cm), ck.bd]), np.log10(ck.ksat_cmh))
    ks_median = 10 ** np.median(np.log10(train.ksat_cmh[train.ksat_cmh > 0]))
    out = []
    for j in np.where(tg[col].isin(held).to_numpy())[0]:
        r = tg.iloc[j]
        theta = st.vg_theta(H, r.thetar, r.thetas, r.alpha_kpa, r.n)
        res = st.estimate(H, theta, ref=ref_bd, n_mc=n_mc, clf=clf_cov,
                          sample_type=r.sample_type, depth=r.depth_cm,
                          bulk_density=r.bd)
        probs = list(res["class_probabilities"].items())
        x = np.array([[np.log10(1 + r.depth_cm), r.bd]])
        k = res["ksat"]
        out.append(dict(
            j=j, design=design, cls=probs[0][0], p1=probs[0][1],
            cls2=probs[1][0] if len(probs) > 1 else None,
            knn=next(iter(res["knn_class_probabilities"])),
            sand=res["fractions"]["sand"], clay=res["fractions"]["clay"],
            ks=k["median_cmh"], ks_lo=k["p5_cmh"], ks_hi=k["p95_cmh"],
            ks_phys=k["physical_matched_cmh"],
            base_cls=base_cls.predict(x)[0], base_ks=10 ** base_ks.predict(x)[0],
            ks_median=ks_median))
    return out


def replicate_floor(g):
    """Spread between replicate Ks cores of one soil (same source, profile
    and depth): the error no predictor of a single measurement can beat."""
    k = g[(g.ksat_cmh > 0) & g.profile_id.notna() & g.depth_cm.notna()]
    rows = []
    for _, x in k.groupby(["source_db", "profile_id", "depth_cm"]):
        if len(x) < 2:
            continue
        v = np.log10(x.ksat_cmh.to_numpy())
        d = np.array([abs(a - b) for a, b in itertools.combinations(v, 2)])
        rows.append((np.median(d), np.mean(d <= LOG2), np.mean(d <= 1)))
    r = np.array(rows)
    pair = np.median(r[:, 0])
    sigma = pair / 0.954          # median |X1 - X2| = 0.954 sigma, lognormal
    from scipy.stats import norm
    return dict(groups=len(r), pair=pair, pair_w2=r[:, 1].mean(), pair_w10=r[:, 2].mean(),
                single=0.674 * sigma, best_w2=2 * norm.cdf(LOG2 / sigma) - 1,
                best_w10=2 * norm.cdf(1 / sigma) - 1)


def metrics(d):
    """Every reported measure for one set of target rows (d has truth columns)."""
    t, p = d.truth.to_numpy(), d.cls.to_numpy()
    m = {"exact": np.mean(p == t),
         "top2": np.mean((p == t) | (d.cls2.to_numpy() == t)),
         "group": np.mean([GROUP[a] == GROUP[b] for a, b in zip(p, t)]),
         "frac_err": np.mean((abs(d.sand - d.t_sand) + abs(d.clay - d.t_clay)
                              + abs((100 - d.sand - d.clay) - (100 - d.t_sand - d.t_clay))) / 3),
         "base_exact": np.mean(d.base_cls.to_numpy() == t),
         "base_group": np.mean([GROUP[a] == GROUP[b] for a, b in zip(d.base_cls, t)])}
    k = d[d.t_ks > 0]
    o = np.log10(k.t_ks.to_numpy())
    for name, col in (("ks", "ks"), ("phys", "ks_phys"), ("base_ks", "base_ks"),
                      ("med_ks", "ks_median")):
        e = np.abs(np.log10(k[col].to_numpy(float)) - o)
        m[f"{name}_w2"] = np.mean(e <= LOG2)
        m[f"{name}_w10"] = np.mean(e <= 1)
        m[f"{name}_factor"] = 10 ** np.median(e)
        m[f"{name}_rho"] = (spearmanr(k[col], k.t_ks).statistic if len(k) > 2 else np.nan)
    m["ks_band"] = np.mean((k.t_ks >= k.ks_lo) & (k.t_ks <= k.ks_hi))
    return m


def boot(d, n=N_BOOT, seed=1):
    """Bootstrap over sources: list of metric dicts."""
    rng = np.random.default_rng(seed)
    srcs = d.source.unique()
    by = {s: d.index[d.source == s] for s in srcs}
    return pd.DataFrame([metrics(d.loc[np.concatenate([by[s] for s in rng.choice(srcs, len(srcs))])])
                         for _ in range(n)])


def ci(b, col):
    lo, hi = np.percentile(b[col].dropna(), [2.5, 97.5])
    return f"[{lo * 100:.0f}–{hi * 100:.0f}]" if col not in ("frac_err",) and "factor" not in col and "rho" not in col \
        else f"[{lo:.2f}–{hi:.2f}]"


def pct(v):
    return f"{v * 100:.0f} %"


def report(d, name, floor):
    m, b = metrics(d), boot(d)
    print(f"\n=== {name}: {len(d)} soils ({int((d.t_ks > 0).sum())} with Ks) from {d.source.nunique()} sources ===")
    print("TEXTURE                       tool (95 % CI)          depth + bulk density alone      chance")
    for key, lab, base, chance in (("exact", "exact class (12)", "base_exact", 1 / 12),
                                   ("top2", "true class in top two", None, 2 / 12),
                                   ("group", "texture group (4)", "base_group", None)):
        diff = b[key] - b[base] if base else None
        bs = (f"{pct(m[base])}  (tool +{(m[key] - m[base]) * 100:.0f} pp, CI "
              f"[{np.percentile(diff, 2.5) * 100:+.0f}, {np.percentile(diff, 97.5) * 100:+.0f}])") if base else "—"
        print(f"  {lab:26s} {pct(m[key]):>5s} {ci(b, key):10s}   {bs:34s} "
              f"{pct(chance) if chance else '25 %*'}")
    print(f"  {'fraction error, points':26s} {m['frac_err']:5.1f} {ci(b, 'frac_err')}")
    print("Ks (measured cores)           within 2x       within 10x      typical factor   rho")
    for key, lab in (("ks", "kNN (the tool's Ks)"), ("phys", "physical (matched)"),
                     ("base_ks", "depth + bulk density alone"), ("med_ks", "reference median (no info)")):
        print(f"  {lab:28s} {pct(m[key + '_w2']):>4s} {ci(b, key + '_w2'):9s} "
              f"{pct(m[key + '_w10']):>4s} {ci(b, key + '_w10'):9s} "
              f"x{m[key + '_factor']:.1f} {ci(b, key + '_factor'):11s} {m[key + '_rho']:.2f}")
    for key in ("ks", "phys"):
        dd = b[f"{key}_w10"] - b["base_ks_w10"]
        d2 = b[f"{key}_rho"] - b["base_ks_rho"]
        print(f"  {key:5s} vs depth + BD alone: within 10x {(m[key + '_w10'] - m['base_ks_w10']) * 100:+.0f} pp "
              f"[{np.percentile(dd, 2.5) * 100:+.0f}, {np.percentile(dd, 97.5) * 100:+.0f}], "
              f"rho {m[key + '_rho'] - m['base_ks_rho']:+.2f} [{np.percentile(d2, 2.5):+.2f}, {np.percentile(d2, 97.5):+.2f}]")
    print(f"  5-95 % band holds {pct(m['ks_band'])} {ci(b, 'ks_band')} (nominal 90 %)")
    print(f"  replicate floor: a single core sits x{10 ** floor['single']:.1f} from its soil's mean; a perfect "
          f"predictor gets ~{pct(floor['best_w2'])} within 2x and ~{pct(floor['best_w10'])} within 10x "
          f"(replicate pairs: {pct(floor['pair_w2'])} / {pct(floor['pair_w10'])}, {floor['groups']} groups)")

    # what one network sees: spread across sources
    per = []
    for s, x in d.groupby("source"):
        if len(x) >= 20:
            mm = metrics(x)
            per.append((s, len(x), mm["exact"], mm["group"],
                        mm["ks_w10"] if (x.t_ks > 0).sum() >= 10 else np.nan))
    per = pd.DataFrame(per, columns=["source", "n", "exact", "group", "ks_w10"])
    q = lambda c: np.nanpercentile(per[c], [10, 50, 90]) * 100
    print(f"Across the {len(per)} sources with >= 20 soils (10th / median / 90th percentile):")
    for c, lab in (("exact", "exact class"), ("group", "texture group"), ("ks_w10", "Ks within 10x")):
        a = q(c)
        print(f"  {lab:14s} {a[0]:3.0f} / {a[1]:3.0f} / {a[2]:3.0f} %")

    # the tool's own confidence signals
    print("Confidence signals:")
    for lo, hi in ((0, .3), (.3, .5), (.5, .7), (.7, 1.01)):
        x = d[(d.p1 >= lo) & (d.p1 < hi)]
        if len(x):
            print(f"  top class probability {lo:.1f}-{min(hi, 1):.1f}: {pct(len(x) / len(d)):>4s} of soils, "
                  f"exact {pct(np.mean(x.cls == x.truth))}, group "
                  f"{pct(np.mean([GROUP[a] == GROUP[b] for a, b in zip(x.cls, x.truth)]))}")
    agree = d.cls == d.knn
    for lab, x in (("classifier and neighbours agree", d[agree]), ("they disagree", d[~agree])):
        print(f"  {lab:32s} {pct(len(x) / len(d)):>4s} of soils, exact {pct(np.mean(x.cls == x.truth))}, "
              f"group {pct(np.mean([GROUP[a] == GROUP[b] for a, b in zip(x.cls, x.truth)]))}")
    k = d[d.t_ks > 0]
    gap = np.abs(np.log10(k.ks / k.ks_phys)) <= 1
    for lab, x in (("kNN and physical Ks within 10x", k[gap]), ("more than 10x apart", k[~gap])):
        e = np.abs(np.log10(x.ks / x.t_ks))
        print(f"  {lab:32s} {pct(len(x) / len(k)):>4s} of soils, Ks within 2x {pct(np.mean(e <= LOG2))}, "
              f"within 10x {pct(np.mean(e <= 1))}")


def frames(n_per_class=150, n_mc=40, reuse=False):
    """Per-target predictions with their truth, by scenario ("source": a new
    network, "profile": an established one), and the Ks replicate floor."""
    # study: the finest unit the metadata support (compilations split into
    # their studies or methods where they can be; verify_network.load)
    from verify_network import load
    g = load()
    pool = g[g.depth_cm.notna() & g.bd.notna()]
    tg = pd.concat([x.sample(min(len(x), n_per_class), random_state=0)
                    for _, x in pool.groupby("texture_class")]).reset_index(drop=True)
    if reuse and os.path.exists(CACHE):
        est = pd.read_csv(CACHE)
    else:
        rng = np.random.default_rng(0)
        profs = np.array(sorted(tg.profile_id.unique()))
        rng.shuffle(profs)
        items = ([("profile", "profile_id", set(x)) for x in np.array_split(profs, 5)]
                 + [("source", "study", {s}) for s in sorted(tg.study.unique())])
        est = pd.DataFrame([r for part in map_folds(fold, items, g, tg, n_mc) for r in part])
        os.makedirs(os.path.dirname(CACHE), exist_ok=True)
        est.to_csv(CACHE, index=False)
    out = {}
    for design in ("source", "profile"):
        d = est[est.design == design].set_index("j").sort_index()
        r = tg.loc[d.index]
        out[design] = d.assign(truth=r.texture_class.to_numpy(), t_sand=r.sand.to_numpy(),
                               t_clay=r.clay.to_numpy(), t_ks=r.ksat_cmh.to_numpy(float),
                               source=r.study.to_numpy()).reset_index(drop=True)
    return out, replicate_floor(g)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    n_per_class = int(args[0]) if args else 150
    n_mc = int(args[1]) if len(args) > 1 else 40
    fr, floor = frames(n_per_class, n_mc, reuse="--reuse" in sys.argv)
    for design, name in (("source", "NEW NETWORK (its study or laboratory never seen), depth and bulk density supplied"),
                         ("profile", "ESTABLISHED NETWORK (its own verified sites in the reference)")):
        report(fr[design], name, floor)
    print("\n* four groups of unequal size; the depth + bulk density model is the fair baseline")


if __name__ == "__main__":
    main()
