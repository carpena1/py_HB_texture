"""How many neighbours should Ks come from, and of which layers?

Until 2026-10 the tool took Ks from the k = 30 nearest reference layers of any
kind, of which only those with a measured Ks voted -- about 14 for a new
source, because half the reference (KSSL) has no Ks. This script scores that
lookup (recomputed here as `legacy_ks`, line for line) against the shipped
one, which searches only layers with a measured Ks, for every k in K_GRID,
and chooses k out of fold, the way ks_physical.AIR_ENTRY_CM was chosen.

Targets as verify_ks_physical.py: per source with >= 40 measured Ks, up to
n_per_source layers stratified by class; curves from the stored vG
parameters, refitted by st.estimate. An undisturbed target passes
sample_type="undisturbed", as a user would. Designs:
  source   the target's whole source hidden          (new site, new source)
  profile  the target's profile hidden               (new site, source in)
  layer    only the target layer hidden              (new depth, other horizons in)
Nested choice: for each scored source, the k whose median |log10 error|,
averaged over the three designs, is lowest on the OTHER sources' targets.

Usage:  python verify_ks_neighbors.py [n_per_source] [n_mc] [--reuse]
Per-target predictions go to figures/cache/ (git-ignored: they include
EU-HYDI layers). Folds run in parallel (verify_common.map_folds).
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

import swcc_texture as st
from verify_common import COLS, map_folds

H = np.arange(0.0, 150.0 + 0.1, 5.0)
MIN_SOURCE = 40
K_GRID = [10, 15, 20, 30, 50, 75, 100, 150, 200, 300, 500]
DESIGNS = ("source", "profile", "layer")
CACHE = os.path.join("figures", "cache", "ks_neighbors.csv")


def legacy_ks(h, theta, ref, n_mc, k=30, seed=0, sample_type=None):
    """Ks as st.estimate computed it before 2026-10: the k nearest layers of
    any kind, those with a measured Ks voting (same fit, draws and weights)."""
    popt, pcov = st.fit_vg(h, theta)
    draws = np.random.default_rng(seed).multivariate_normal(popt, pcov, size=n_mc)
    draws = np.clip(draws, [0.0, 0.05, -4.0, np.log10(0.01)],
                    [0.5, 1.0, 1.5, np.log10(10.0)])
    ks_type = "undisturbed" if sample_type == "undisturbed" else None
    vals, wts = [], []
    for tr, ts, la, ln1 in draws:
        idx, w = ref.neighbors(tr, ts, 10.0 ** la, 1.0 + 10.0 ** ln1, k,
                               sample_type=ks_type)
        w = w * ref.class_weight[idx] * ref.row_weight[idx]
        w = w / w.sum()
        ok = ref.has_ksat[idx]
        vals.extend(ref.ksat[idx][ok]); wts.extend(w[ok])
    if not vals:
        return np.nan, np.nan, np.nan
    lk, wts = np.log10(vals), np.array(wts)
    return tuple(10.0 ** st._wpercentile(lk, wts, q) for q in (50, 5, 95))


def run_source(s, df, n_per_source, n_mc):
    ref = st.GshpReference(df=df[COLS + ["source_db", "sample_type"]].reset_index(drop=True))
    src, prof, lid = df.source_db.to_numpy(), df.profile_id.to_numpy(), df.layer_id.to_numpy()
    sub = df[(src == s) & ref.has_ksat]
    tg = pd.concat([g.sample(min(len(g), max(1, n_per_source // 12)), random_state=0)
                    for _, g in sub.groupby("texture_class")])
    out = []
    for row in tg.itertuples():
        theta = st.vg_theta(H, row.thetar, row.thetas, row.alpha_kpa, row.n)
        stype = "undisturbed" if row.sample_type == "undisturbed" else None
        for design, blocked in zip(DESIGNS, (src == s, prof == row.profile_id,
                                             lid == row.layer_id)):
            ref.set_excluded(blocked)
            rec = dict(source=s, design=design, layer_id=row.layer_id, obs=row.ksat_cmh)
            rec["legacy"], rec["legacy_lo"], rec["legacy_hi"] = legacy_ks(
                H, theta, ref, n_mc, sample_type=stype)
            for kk in K_GRID:
                k = st.estimate(H, theta, ref=ref, n_mc=n_mc, sample_type=stype, k_ks=kk)["ksat"]
                rec.update({f"k{kk}": k["median_cmh"], f"k{kk}_lo": k["p5_cmh"],
                            f"k{kk}_hi": k["p95_cmh"]})
            out.append(rec)
    return pd.DataFrame(out)


def ae(pred, obs):
    return np.abs(np.log10(np.asarray(pred, float)) - np.log10(np.asarray(obs, float)))


def scores(label, x, col):
    p = x[col].to_numpy(float)
    o = x.obs.to_numpy(float)
    e = np.log10(p) - np.log10(o)
    band = np.mean((o >= x[f"{col}_lo"]) & (o <= x[f"{col}_hi"]))
    print(f"    {label:26s} median error {np.median(abs(e)):.3f} dex  bias {np.median(e):+.2f}"
          f"  within 2x {100 * np.mean(abs(e) <= np.log10(2)):3.0f}%"
          f"  10x {100 * np.mean(abs(e) <= 1):3.0f}%  band {100 * band:3.0f}%"
          f"  rho {spearmanr(p, o).statistic:.3f}")
    return abs(e)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    n_per_source = int(args[0]) if args else 200
    n_mc = int(args[1]) if len(args) > 1 else 20
    if "--reuse" in sys.argv and os.path.exists(CACHE):
        d = pd.read_csv(CACHE)
    else:
        df = st.load_reference_df(st.DEFAULT_REFERENCE).reset_index(drop=True)
        df["source_db"] = df.source_db.fillna("KSSL")
        df["profile_id"] = df.profile_id.fillna("solo_" + df.layer_id.astype(str))
        has = df.ksat_cmh.notna() & (df.ksat_cmh > 0)
        sizes = df[has].source_db.value_counts()
        sources = [s for s in sizes.index if sizes[s] >= MIN_SOURCE]
        print(f"reference {len(df)} layers, {has.sum()} with a measured Ks; "
              f"{len(sources)} sources scored, n_mc={n_mc}")
        d = pd.concat(map_folds(run_source, sources, df, n_per_source, n_mc),
                      ignore_index=True)
        os.makedirs(os.path.dirname(CACHE), exist_ok=True)
        d.to_csv(CACHE, index=False)

    def crit(x, kk):
        return np.mean([np.median(ae(x[x.design == g][f"k{kk}"], x[x.design == g].obs))
                        for g in DESIGNS])

    print("\nmedian |error| by k, Ks-bearing layers only (pooled over sources):")
    for g in DESIGNS:
        x = d[d.design == g]
        print(f"  {g:8s} legacy {np.median(ae(x.legacy, x.obs)):.3f} | "
              + "  ".join(f"k{kk} {np.median(ae(x[f'k{kk}'], x.obs)):.3f}" for kk in K_GRID))
    chosen = {s: min(K_GRID, key=lambda kk: crit(d[d.source != s], kk))
              for s in d.source.unique()}
    print("\nk chosen out of fold (mean over the three designs): "
          + ", ".join(f"{k} {v}x" for k, v in pd.Series(chosen).value_counts().sort_index().items())
          + f"; the tool uses KS_K = {st.KS_K}")

    for g in DESIGNS:
        x = d[d.design == g].reset_index(drop=True)
        for suffix in ("", "_lo", "_hi"):
            x[f"nested{suffix}"] = [r[f"k{chosen[r.source]}{suffix}"] for _, r in x.iterrows()]
        print(f"\n=== {g} design: {len(x)} targets, {x.source.nunique()} sources ===")
        e_old = scores("legacy (30, any layer)", x, "legacy")
        scores("Ks-bearing, k=30", x, "k30")
        e_nest = scores("Ks-bearing, k out of fold", x, "nested")
        scores(f"Ks-bearing, k={st.KS_K} (tool)", x, f"k{st.KS_K}")
        dd = e_nest - e_old
        per = x.assign(a=e_old, b=e_nest).groupby("source")[["a", "b"]].median()
        print(f"    out of fold vs legacy: mean {dd.mean():+.3f} dex, better {np.mean(dd < 0) * 100:.0f}% /"
              f" worse {np.mean(dd > 0) * 100:.0f}%, Wilcoxon p={wilcoxon(dd[dd != 0]).pvalue:.2g};"
              f" sources better {(per.b < per.a).sum()} of {len(per)}")


if __name__ == "__main__":
    main()
