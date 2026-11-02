"""Does freeing m from n improve texture and Ks prediction?

The Mualem constraint m = 1 - 1/n ties the two shape parameters together, so
the reference is indexed on four numbers of which two are redundant by
construction (corr(log10(n-1), 1-1/n) = 0.993). Freeing m gives the dry limb
of the curve -- the part most closely tied to texture -- its own coordinate;
in the rebuilt GSHP table the correlation drops to -0.38, so it is a genuinely
new axis rather than a restatement of n.

DESIGN. Each target's synthetic curve is generated from its UNCONSTRAINED
parameters, because those are the closer description of the measured points
the row was fitted to. Both arms then see that identical curve and refit it
their own way:

    A  4 features, Mualem refit, matched against the Mualem reference
    B  5 features, free-m refit, matched against the free-m reference

Generating from the Mualem parameters instead would hand arm B a curve that is
Mualem-shaped by construction, which would rig the comparison against it.

The evaluation grid spans 0.1 to 1500 kPa, not the 0-150 kPa grid the older
benchmarks use: m governs the dry limb, and a grid stopping at 150 kPa cannot
see it.

Both objectives are reported (texture and Ksat) for both models (kNN and the
hybrid GBM). Profile-level leave-one-out.

Usage:  python verify_free_m.py [n_per_class] [n_mc]
"""

import sys

import numpy as np
import pandas as pd

import ksat_metrics as km
import swcc_texture as st
from verify_by_class import mcnemar
from verify_groups import GROUP, GROUP_ORDER
from verify_tau import prf
from verify_variants import ORDER

# Wide grid: m acts on the dry end, so the curve must reach it.
H = np.concatenate([[0.0], np.logspace(-1, np.log10(1500.0), 25)])


N_FOLDS = 5


def run(df, targets, use_m, n_mc, hybrid, folds):
    """Profile-grouped k-fold. Both halves of the hybrid must be blind to the
    target: the kNN reference and the GBM's training set are the SAME reduced
    frame. Training the classifier once on the full table instead lets it
    memorise every target, which inflates the hybrid arms and flatters the
    five-feature one most, since five coordinates identify a row more
    precisely than four."""
    cls = np.empty(len(targets), dtype=object)
    med = np.full(len(targets), np.nan)
    lo = np.full(len(targets), np.nan)
    hi = np.full(len(targets), np.nan)
    tp = targets.profile_id.to_numpy()
    for held in folds:
        train = df[~df.profile_id.isin(held)].reset_index(drop=True)
        ref = st.GshpReference(df=train, use_m=use_m)
        clf = st.TextureGBM(df=train, use_m=use_m) if hybrid else None
        sel = np.where(np.isin(tp, list(held)))[0]
        for j in sel:
            row = targets.iloc[j]
            # curve generated from the unconstrained parameters (docstring)
            theta = st.vg_theta(H, row.thetar_m, row.thetas_m,
                                row.alpha_kpa_m, row.n_m, row.m)
            r = st.estimate(H, theta, ref=ref, n_mc=n_mc, clf=clf)
            cls[j] = r["texture_class"]
            k = r["ksat"]
            med[j], lo[j], hi[j] = (k["median_cmh"], k["p5_cmh"], k["p95_cmh"])
    return np.array(cls), med, lo, hi


def report(truth, tgroup, obs, res):
    print(f"\n{'variant':<22s} {'overall':>8s} {'macroR':>7s} {'macroP':>7s} "
          f"{'macroF1':>8s} {'group':>7s}")
    print("-" * 62)
    for name, out in res.items():
        p = out[0]
        R, P, F = prf(truth, p)
        grp = np.mean([GROUP[a] == GROUP[t] for a, t in zip(p, truth)])
        print(f"{name:<22s} {(p == truth).mean()*100:7.1f}% {R*100:6.1f}% "
              f"{P*100:6.1f}% {F*100:7.1f}% {grp*100:6.1f}%")

    keys = list(res)
    for i in range(0, len(keys), 2):
        a, b = keys[i], keys[i + 1]
        ca, cb = res[a][0] == truth, res[b][0] == truth
        ga = np.array([GROUP[x] == GROUP[t] for x, t in zip(res[a][0], truth)])
        gb = np.array([GROUP[x] == GROUP[t] for x, t in zip(res[b][0], truth)])
        print(f"\n  {b} vs {a}:  class {(cb.mean()-ca.mean())*100:+.1f}pp "
              f"p={mcnemar(ca, cb):.4f}   group "
              f"{(gb.mean()-ga.mean())*100:+.1f}pp p={mcnemar(ga, gb):.4f}")

    print("\nPER-CLASS RECALL")
    hdr = "".join(f"{k:>12s}" for k in res)
    print(f"{'class':<18s} {'n':>5s}{hdr}")
    for cl in ORDER:
        m = truth == cl
        if not m.any():
            continue
        vals = "".join(f"{(res[k][0][m] == cl).mean()*100:11.1f}%" for k in res)
        print(f"{cl:<18s} {m.sum():5d}{vals}")

    print(f"\n{'group':<18s} {'n':>5s}{hdr}")
    for g in GROUP_ORDER:
        m = tgroup == g
        if not m.any():
            continue
        vals = "".join(
            f"{np.mean([GROUP[x] == GROUP[t] for x, t in zip(res[k][0][m], truth[m])])*100:11.1f}%"
            for k in res)
        print(f"{g:<18s} {m.sum():5d}{vals}")

    km.report("KSAT", [(k, km.ksat_scores(obs, v[1], v[2], v[3]))
                       for k, v in res.items()])

    # Paired test on Ksat: |log10 error| per soil, Wilcoxon signed-rank on the
    # soils both arms predicted. Ks always comes from the kNN, so the two
    # hybrid arms reproduce the two kNN arms exactly here.
    from scipy.stats import wilcoxon
    keys = list(res)
    print("\nPAIRED KSAT (|log10 error| per soil, Wilcoxon signed-rank)")
    for i in range(0, len(keys), 2):
        a, b = keys[i], keys[i + 1]
        ea = np.abs(np.log10(res[a][1]) - np.log10(obs))
        eb = np.abs(np.log10(res[b][1]) - np.log10(obs))
        ok = np.isfinite(ea) & np.isfinite(eb)
        if ok.sum() < 10:
            continue
        stat, pv = wilcoxon(ea[ok], eb[ok])
        print(f"  {b} vs {a}: n={ok.sum()}  median |err| "
              f"{np.median(ea[ok]):.3f} -> {np.median(eb[ok]):.3f} dex  "
              f"mean {ea[ok].mean():.3f} -> {eb[ok].mean():.3f}  p={pv:.4f}")


def main():
    n_per_class = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    n_mc = int(sys.argv[2]) if len(sys.argv) > 2 else 50

    df = st.load_reference_df()
    df["profile_id"] = df.profile_id.fillna("solo_" + df.layer_id.astype(str))
    pool = df[df[["thetar_m", "thetas_m", "alpha_kpa_m", "n_m", "m"]]
              .notna().all(axis=1)]
    targets = pd.concat([g.sample(min(len(g), n_per_class), random_state=0)
                         for _, g in pool.groupby("texture_class")])
    truth = targets.texture_class.to_numpy()
    tgroup = np.array([GROUP[t] for t in truth])
    obs = targets.ksat_cmh.to_numpy(float)

    rng = np.random.default_rng(0)
    profs = targets.profile_id.unique().copy()
    rng.shuffle(profs)
    folds = [set(x) for x in np.array_split(profs, N_FOLDS)]

    print(f"reference {len(df)} layers (merged: GSHP + KSSL + Hohenbrink)")
    print(f"targets {len(targets)}; profile-grouped {N_FOLDS}-fold "
          f"(kNN and GBM share each fold's training frame); n_mc={n_mc}; "
          f"grid {H[1]:.2g}-{H[-1]:.0f} kPa")
    print(f"free m: median {df.m.median():.3f} vs Mualem "
          f"{(1 - 1 / df.n).median():.3f}")

    res = {"kNN 4 (Mualem)": run(df, targets, False, n_mc, False, folds),
           "kNN 5 (+m)": run(df, targets, True, n_mc, False, folds),
           "hybrid 4 (Mualem)": run(df, targets, False, n_mc, True, folds),
           "hybrid 5 (+m)": run(df, targets, True, n_mc, True, folds)}
    report(truth, tgroup, obs, res)


if __name__ == "__main__":
    main()
