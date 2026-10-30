"""Do the soils whose texture may be on a < 1 mm sieving basis help or hurt the
rest of the reference?

Sources in swcc_texture.SIEVE_1MM_SOURCES may give sand as 50-1000 um, not
50-2000 um (Polish laboratories before about 2008; the Russian standard,
per A. Nemes, 2026-10). The project's standard targets -- reference layers
NOT in those sources, 150 per class -- are scored twice on identical folds,
with the flagged layers in the reference (the default) and taken out, under
both hold-out designs (new site: profile held out; new source: the target's
whole contributing source held out). Ks is scored from the kNN alone with
each target's source held out, as in verify_join_effect.py. The predicted
sand fraction is kept too, since a < 1 mm basis would mainly shift sand.

Usage:  python verify_sieve_1mm.py [n_per_class] [n_mc]
"""

import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import swcc_texture as st
from verify_common import map_folds, mcnemar
from verify_groups import GROUP
from verify_join_effect import H, N_FOLDS, REF_COLS, curve, prep


def score(job, ref_df, tg, n_mc):
    """One reference (and classifier) per fold; (class, sand %, Ks cm/h)."""
    col, folds, clf = job
    pred = np.empty(len(tg), dtype=object)
    sand = np.full(len(tg), np.nan)
    ks = np.full(len(tg), np.nan)
    for held in folds:
        sel = np.where(tg[col].isin(held).to_numpy())[0]
        if not len(sel):
            continue
        train = ref_df[~ref_df[col].isin(held)]
        ref = st.GshpReference(df=train[REF_COLS].reset_index(drop=True))
        model = st.TextureGBM(df=train.reset_index(drop=True)) if clf else None
        for j in sel:
            r = tg.iloc[j]
            out = st.estimate(H, curve(r), ref=ref, n_mc=n_mc, clf=model,
                              sample_type=r.sample_type)
            pred[j] = out["texture_class"]
            sand[j] = out["fractions"]["sand"]
            ks[j] = out["ksat"]["median_cmh"]
    return pred, sand, ks


def run(job, refs, targets, n_mc):
    arm, kind, col, folds, clf = job
    return score((col, folds, clf), refs[arm], targets[kind], n_mc)


def main():
    n_per_class = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    n_mc = int(sys.argv[2]) if len(sys.argv) > 2 else 30

    both = prep(st.load_reference_df())
    flag = both.source_db.isin(st.SIEVE_1MM_SOURCES)
    base = both[~flag].reset_index(drop=True)
    print(f"with the flagged sources {len(both)} layers, without {len(base)} "
          f"(-{int(flag.sum())}: {both[flag].source_db.value_counts().to_dict()})",
          flush=True)

    tg = pd.concat([g.sample(min(len(g), n_per_class), random_state=0)
                    for _, g in base.groupby("texture_class")]).reset_index(drop=True)
    pool = base[base.ksat_cmh.notna() & (base.ksat_cmh > 0)]
    tk = pd.concat([g.sample(min(len(g), n_per_class), random_state=0)
                    for _, g in pool.groupby("texture_class")]).reset_index(drop=True)
    rng = np.random.default_rng(0)
    profs = np.array(sorted(tg.profile_id.unique())); rng.shuffle(profs)
    designs = {"new site": ("profile_id", [set(x) for x in np.array_split(profs, N_FOLDS)]),
               "new source": ("source_db", [{s} for s in sorted(tg.source_db.unique())])}
    jobs = [(arm, "tex", col, folds, True) for col, folds in designs.values()
            for arm in ("with", "without")]
    jobs += [(arm, "ks", "source_db", [{s} for s in sorted(tk.source_db.unique())], False)
             for arm in ("with", "without")]
    print(f"targets {len(tg)} ({n_per_class}/class), Ks targets {len(tk)}, n_mc={n_mc}\n",
          flush=True)
    out = map_folds(run, jobs, {"with": both, "without": base}, {"tex": tg, "ks": tk},
                    n_mc, threads=3)
    res = {(j[0], j[2] if j[1] == "tex" else "ks"): o for j, o in zip(jobs, out)}

    truth = tg.texture_class.to_numpy()
    grp = lambda p: np.array([GROUP[a] == GROUP[t] for a, t in zip(p, truth)])
    print("EFFECT OF KEEPING THE FLAGGED SOURCES (with minus without)")
    for design, (col, _) in designs.items():
        (pw, sw, _), (po, so, _) = res[("with", col)], res[("without", col)]
        a, b = po == truth, pw == truth
        print(f"  {design:<11s} exact {a.mean()*100:5.2f} -> {b.mean()*100:5.2f} % "
              f"({(b.mean()-a.mean())*100:+.2f} pp, p={mcnemar(a, b):.3f})   group "
              f"{(grp(pw).mean()-grp(po).mean())*100:+.2f} pp (p={mcnemar(grp(po), grp(pw)):.3f})"
              f"   answers changed {int((pw != po).sum())}/{len(truth)}")
        e_w, e_o = sw - tg.sand.to_numpy(), so - tg.sand.to_numpy()
        print(f"  {'':<11s} sand bias {np.mean(e_o):+.2f} -> {np.mean(e_w):+.2f} points, "
              f"mean |error| {np.mean(np.abs(e_o)):.2f} -> {np.mean(np.abs(e_w)):.2f} "
              f"(Wilcoxon p={wilcoxon(np.abs(e_w), np.abs(e_o)).pvalue:.3f})")
        for g in ("Sandy", "Loamy", "Silty", "Clayey"):
            m = np.array([GROUP[t] == g for t in truth])
            print(f"      {g:<7s} n={m.sum():4d}  exact {a[m].mean()*100:5.1f} -> "
                  f"{b[m].mean()*100:5.1f} %   sand bias {np.mean(e_o[m]):+5.1f} -> "
                  f"{np.mean(e_w[m]):+5.1f}")

    obs = np.log10(tk.ksat_cmh.to_numpy(float))
    err = {arm: np.abs(np.log10(res[(arm, "ks")][2]) - obs) for arm in ("with", "without")}
    ok = np.isfinite(err["with"]) & np.isfinite(err["without"])
    d = err["with"][ok] - err["without"][ok]
    print(f"\nKs ({int(ok.sum())} targets): median |log10 err| "
          f"{np.median(err['without'][ok]):.4f} -> {np.median(err['with'][ok]):.4f} dex, "
          f"within 10x {np.mean(err['without'][ok] < 1)*100:.1f} -> "
          f"{np.mean(err['with'][ok] < 1)*100:.1f} %, changed {int((d != 0).sum())}"
          + (f", Wilcoxon p={wilcoxon(d[d != 0]).pvalue:.3f}" if (d != 0).any() else ""))


if __name__ == "__main__":
    main()
