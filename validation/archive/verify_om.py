"""Does organic carbon earn the reference rows it costs?

Organic matter was dropped early on because GSHP carries it for only ~15 % of
layers. KSSL carries it for ~99 %, so the merged reference is now 32 %
covered (4,023 of 12,526) and the question is worth reopening.

The catch is that using it as a kNN feature forces every reference row without
it to be discarded. So there are two separate effects, and a naive test
conflates them:

  A  OM-only reference, 4 features  -- the cost of shrinking, isolated
  B  OM-only reference, 5 features  -- A plus the covariate
  C  full merged reference, 4 features -- what the tool does today

B vs A isolates whether the covariate carries information.
B vs C is the decision that actually matters: is the whole package better
than today's default?

Targets are soils that themselves have organic carbon, since a user could not
supply the covariate otherwise. Profile-level leave-one-out; kNN.
Texture and Ksat are both reported.

Usage:  python verify_om.py [n_per_class] [n_mc]
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

H = np.arange(0.0, 150.0 + 0.1, 5.0)
COLS = ["layer_id", "profile_id", "texture_class", "alpha_kpa", "n", "thetar",
        "thetas", "sand", "silt", "clay", "ksat_cmh", "depth_cm", "oc"]


def run(ref_df, targets, use_om, n_mc):
    ref = st.GshpReference(df=ref_df.reset_index(drop=True), use_om=use_om)
    cls, med, lo, hi = [], [], [], []
    for row in targets.itertuples():
        theta = st.vg_theta(H, row.thetar, row.thetas, row.alpha_kpa, row.n)
        ref.set_excluded_profile(row.profile_id)
        r = st.estimate(H, theta, ref=ref, n_mc=n_mc,
                        om=row.oc if use_om else None)
        cls.append(r["texture_class"])
        k = r["ksat"]
        med.append(k["median_cmh"]); lo.append(k["p5_cmh"]); hi.append(k["p95_cmh"])
    return (np.array(cls), np.array(med, float), np.array(lo, float),
            np.array(hi, float))


def main():
    n_per_class = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    n_mc = int(sys.argv[2]) if len(sys.argv) > 2 else 50

    full = st.load_reference_df("merged")[COLS]
    om_only = full[full.oc.notna()]
    pool = om_only[om_only.profile_id.notna()]
    targets = pd.concat([g.sample(min(len(g), n_per_class), random_state=0)
                         for _, g in pool.groupby("texture_class")])
    truth = targets.texture_class.to_numpy()
    tgroup = np.array([GROUP[t] for t in truth])
    obs = targets.ksat_cmh.to_numpy(float)

    print(f"merged reference {len(full)} layers; with organic carbon "
          f"{len(om_only)} ({len(om_only)/len(full)*100:.1f}%)")
    print(f"targets {len(targets)} (all carry organic carbon); "
          f"profile-level LOO; kNN; n_mc={n_mc}")
    print(f"target organic carbon: median {targets.oc.median():.2f}%, "
          f"p90 {targets.oc.quantile(.9):.2f}%\n")

    res = {
        "A OM-ref, no OM": run(om_only, targets, False, n_mc),
        "B OM-ref, +OM": run(om_only, targets, True, n_mc),
        "C full ref (today)": run(full, targets, False, n_mc),
    }

    print(f"{'variant':<20s} {'overall':>8s} {'macroR':>7s} {'macroP':>7s} "
          f"{'macroF1':>8s} {'group':>7s} {'ref rows':>9s}")
    print("-" * 68)
    for name, out in res.items():
        p = out[0]
        R, P, F = prf(truth, p)
        grp = np.mean([GROUP[a] == GROUP[t] for a, t in zip(p, truth)])
        rows = len(om_only) if name.startswith(("A", "B")) else len(full)
        print(f"{name:<20s} {(p == truth).mean()*100:7.1f}% {R*100:6.1f}% "
              f"{P*100:6.1f}% {F*100:7.1f}% {grp*100:6.1f}% {rows:9d}")

    a, b, c = (res[k][0] == truth for k in res)
    print(f"\nB vs A  (does the covariate carry information?)  "
          f"{(b.mean()-a.mean())*100:+.1f}pp   p = {mcnemar(a, b):.4f}")
    print(f"B vs C  (is the whole package better than today?) "
          f"{(b.mean()-c.mean())*100:+.1f}pp   p = {mcnemar(c, b):.4f}")
    print(f"A vs C  (cost of shrinking the reference)         "
          f"{(a.mean()-c.mean())*100:+.1f}pp   p = {mcnemar(c, a):.4f}")

    print("\nPER-CLASS RECALL")
    print(f"{'class':<18s} {'n':>5s} {'A':>7s} {'B':>7s} {'C':>7s} "
          f"{'B-A':>7s} {'B-C':>7s}")
    for cl in ORDER:
        m = truth == cl
        if not m.any():
            continue
        va, vb, vc = ((res[k][0][m] == cl).mean() for k in res)
        print(f"{cl:<18s} {m.sum():5d} {va*100:6.1f}% {vb*100:6.1f}% "
              f"{vc*100:6.1f}% {(vb-va)*100:+6.1f}pp {(vb-vc)*100:+6.1f}pp")

    print(f"\n{'group':<18s} {'n':>5s} {'A':>7s} {'B':>7s} {'C':>7s} "
          f"{'B-A':>7s} {'B-C':>7s}")
    for g in GROUP_ORDER:
        m = tgroup == g
        if not m.any():
            continue
        va, vb, vc = (np.mean([GROUP[x] == GROUP[t]
                               for x, t in zip(res[k][0][m], truth[m])])
                      for k in res)
        print(f"{g:<18s} {m.sum():5d} {va*100:6.1f}% {vb*100:6.1f}% "
              f"{vc*100:6.1f}% {(vb-va)*100:+6.1f}pp {(vb-vc)*100:+6.1f}pp")

    km.report("KSAT", [(k, km.ksat_scores(obs, v[1], v[2], v[3]))
                       for k, v in res.items()])


if __name__ == "__main__":
    main()
