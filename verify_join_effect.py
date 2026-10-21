"""Does adding a candidate table make the predictions for the REST worse?

The admission test every external set has had to pass before joining the
default reference (README, "Andalusian soils"). The project's standard
targets -- reference layers that are NOT in the candidate, 150 per class --
are scored twice on identical folds, with the candidate out of the reference
and in it, under both hold-out designs:

  new site    the target's profile is held out of the reference and of the
              classifier's training data (profile-grouped folds)
  new source  the target's whole contributing laboratory is held out of both

Ks is scored separately, from the kNN alone (no classifier), with each
target's source held out. A candidate with no measured Ks of its own can
still move Ks slightly: the neighbour vote weights by inverse class
frequency, so extra rows of one class reweight the donors that are selected.

Run it BEFORE registering the candidate in REFERENCE_SETS -- it asserts the
candidate is absent from the default, because that is its baseline.

Result for Lebrija (2026-10, 1,733 targets / 1,660 with Ks, n_mc=30):
no harm, nothing significant. Exact class +0.87 pp for a new site
(p=0.176) and +0.17 pp for a new source (p=0.839); texture group +1.04
(p=0.092) and -0.29 (p=0.694); Ks median |log10 error| +0.0009 dex with
107 of 1,660 predictions changed at all (Wilcoxon p=0.872). Admitted.

Usage:  python verify_join_effect.py <candidate> [n_per_class] [n_mc]
          <candidate>  a name in swcc_texture.REFERENCE_SETS, e.g. lebrija
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import swcc_texture as st
from verify_common import mcnemar, COLS
from verify_groups import GROUP

H = np.concatenate([[0.0], np.logspace(-1, np.log10(1500.0), 25)])
N_FOLDS = 5
REF_COLS = COLS + ["source_db", "sample_type"]


def curve(r):
    return st.vg_theta(H, r.thetar, r.thetas, r.alpha_kpa, r.n)


def prep(d):
    d = d.copy()
    d["source_db"] = d.source_db.fillna("KSSL")
    d["profile_id"] = d.profile_id.fillna("solo_" + d.layer_id.astype(str))
    d["sample_type"] = d.sample_type.fillna("unknown")
    return d


def score(ref_df, tg, folds, col, n_mc, clf=True):
    """One reference (and classifier) per fold, not per target."""
    pred = np.empty(len(tg), dtype=object)
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
            ks[j] = out["ksat"]["median_cmh"]
    return pred, ks


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in st.REFERENCE_SETS:
        raise SystemExit("usage: verify_join_effect.py <candidate> "
                         "[n_per_class] [n_mc]\n  candidate in "
                         f"{sorted(st.REFERENCE_SETS)}")
    name = sys.argv[1]
    n_per_class = int(sys.argv[2]) if len(sys.argv) > 2 else 150
    n_mc = int(sys.argv[3]) if len(sys.argv) > 3 else 30

    base = prep(st.load_reference_df())
    cand = prep(st.load_reference_df(name))
    src = set(cand.source_db.dropna())
    assert not base.source_db.isin(src).any(), (
        f"{name} is already in the default reference -- run this before "
        "registering it in REFERENCE_SETS")
    both = pd.concat([base, cand], ignore_index=True)
    print(f"without {name} {len(base)} layers, with it {len(both)} "
          f"(+{len(cand)})", flush=True)

    tg = pd.concat([g.sample(min(len(g), n_per_class), random_state=0)
                    for _, g in base.groupby("texture_class")]
                   ).reset_index(drop=True)
    truth = tg.texture_class.to_numpy()
    rng = np.random.default_rng(0)
    profs = np.array(sorted(tg.profile_id.unique())); rng.shuffle(profs)
    designs = {
        "new site": ("profile_id",
                     [set(x) for x in np.array_split(profs, N_FOLDS)]),
        "new source": ("source_db", [{s} for s in sorted(tg.source_db.unique())]),
    }
    print(f"targets {len(tg)} ({n_per_class}/class), n_mc={n_mc}; "
          f"{len(designs['new source'][1])} source folds\n", flush=True)

    res = {}
    for design, (col, folds) in designs.items():
        for arm, ref_df in (("without", base), ("with", both)):
            res[(arm, design)] = score(ref_df, tg, folds, col, n_mc)
            p = res[(arm, design)][0]
            g = np.mean([GROUP[a] == GROUP[b] for a, b in zip(p, truth)])
            print(f"  {arm:<8s} {design:<11s} exact "
                  f"{np.mean(p == truth)*100:5.2f} % group {g*100:5.2f} %",
                  flush=True)

    print("\nEFFECT ON THE REST OF THE REFERENCE (with minus without)")
    for design in designs:
        pa, pb = res[("without", design)][0], res[("with", design)][0]
        a, b = pa == truth, pb == truth
        ga = np.array([GROUP[x] == GROUP[t] for x, t in zip(pa, truth)])
        gb = np.array([GROUP[x] == GROUP[t] for x, t in zip(pb, truth)])
        print(f"  {design:<11s} exact {(b.mean()-a.mean())*100:+5.2f} pp "
              f"(p={mcnemar(a, b):.3f})   group "
              f"{(gb.mean()-ga.mean())*100:+5.2f} pp "
              f"(p={mcnemar(ga, gb):.3f})   answers changed "
              f"{int((pa != pb).sum())}/{len(truth)}")

    # Ks: kNN only, each target's source held out.
    pool = base[base.ksat_cmh.notna() & (base.ksat_cmh > 0)]
    tk = pd.concat([g.sample(min(len(g), n_per_class), random_state=0)
                    for _, g in pool.groupby("texture_class")]
                   ).reset_index(drop=True)
    obs = np.log10(tk.ksat_cmh.to_numpy(float))
    print(f"\nKs: {len(tk)} targets with a measured Ks", flush=True)
    err = {}
    for arm, ref_df in (("without", base), ("with", both)):
        folds = [{s} for s in sorted(tk.source_db.unique())]
        pred = score(ref_df, tk, folds, "source_db", n_mc, clf=False)[1]
        err[arm] = np.abs(np.log10(pred) - obs)
        ok = np.isfinite(err[arm])
        print(f"  {arm:<8s} median |log10 err| {np.median(err[arm][ok]):.4f} dex"
              f"  within 2x {np.mean(err[arm][ok] < np.log10(2))*100:.1f} %"
              f"  within 10x {np.mean(err[arm][ok] < 1)*100:.1f} %", flush=True)
    ok = np.isfinite(err["without"]) & np.isfinite(err["with"])
    d = err["with"][ok] - err["without"][ok]
    print(f"  with minus without: median "
          f"{np.median(err['with'][ok]) - np.median(err['without'][ok]):+.4f} dex,"
          f" changed {int((d != 0).sum())}/{int(ok.sum())}"
          + (f", Wilcoxon p={wilcoxon(d[d != 0]).pvalue:.3f}"
             if (d != 0).any() else ""))


if __name__ == "__main__":
    main()
