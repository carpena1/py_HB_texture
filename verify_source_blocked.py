"""How much of the reported accuracy survives when the tool meets a soil from
a laboratory it has never seen?

Every accuracy figure in this project so far comes from leave-one-out with the
target's *profile* removed. That rules out site leakage, but not source
leakage: GSHP is a compilation, and one contributor (Florida_database) supplies
57 % of the curated layers. Under profile-level LOO a Florida target is still
matched against thousands of other Florida layers, measured by the same lab
with the same protocol on the same regional soils. Real users are not in that
position.

This script holds out an entire contributing database at once and predicts its
layers from the remaining sources. The same targets are also predicted under
ordinary profile-level LOO, so the two are paired and the gap is the leakage.

Both objectives are reported: texture (exact class and aggregated group) and
Ksat (ksat_metrics).

Usage:  python verify_source_blocked.py [n_per_source] [n_mc]
"""

import sys

import numpy as np
import pandas as pd

import ksat_metrics as km
import swcc_texture as st
from verify_by_class import mcnemar
from verify_groups import GROUP
from verify_variants import COLS

MIN_SOURCE = 80          # skip contributors too small to score meaningfully
H = np.arange(0.0, 150.0 + 0.1, 5.0)
REF_COLS = COLS + ["source_db"]


def predict(ref, targets, blocked_mask, n_mc):
    """Predict each target. If blocked_mask is given, the whole source is
    hidden; otherwise only the target's own profile is."""
    cls, med, lo, hi = [], [], [], []
    for row in targets.itertuples():
        theta = st.vg_theta(H, row.thetar, row.thetas, row.alpha_kpa, row.n)
        if blocked_mask is None:
            ref.set_excluded_profile(row.profile_id)
        else:
            ref.set_excluded(blocked_mask)
        r = st.estimate(H, theta, ref=ref, n_mc=n_mc)
        cls.append(r["texture_class"])
        k = r["ksat"]
        med.append(k["median_cmh"]); lo.append(k["p5_cmh"]); hi.append(k["p95_cmh"])
    return (np.array(cls), np.array(med, float), np.array(lo, float),
            np.array(hi, float))


def main():
    n_per_source = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    n_mc = int(sys.argv[2]) if len(sys.argv) > 2 else 60

    gshp = pd.read_csv(st.REFERENCE_CSV)
    ref = st.GshpReference(df=gshp[REF_COLS].reset_index(drop=True))
    src_all = gshp.source_db.to_numpy()

    sizes = gshp.source_db.value_counts()
    sources = [s for s in sizes.index if sizes[s] >= MIN_SOURCE]
    print(f"reference {len(gshp)} layers from {len(sizes)} contributing "
          f"databases; scoring the {len(sources)} with >= {MIN_SOURCE} layers")
    print(f"n_mc={n_mc}, up to {n_per_source} targets per source\n")

    rows, pooled = [], {"prof": [], "block": [], "gprof": [], "gblock": [],
                        "obs": [], "kb": [[], [], []]}
    for s in sources:
        sub = gshp[gshp.source_db == s]
        # stratify within the source so its own class mix does not dominate
        tg = pd.concat([g.sample(min(len(g), max(1, n_per_source // 12)),
                                 random_state=0)
                        for _, g in sub.groupby("texture_class")])
        truth = tg.texture_class.to_numpy()
        mask = (src_all == s)

        p_prof = predict(ref, tg, None, n_mc)
        p_block = predict(ref, tg, mask, n_mc)
        a, b = p_prof[0] == truth, p_block[0] == truth
        ga = np.array([GROUP[p] == GROUP[t] for p, t in zip(p_prof[0], truth)])
        gb = np.array([GROUP[p] == GROUP[t] for p, t in zip(p_block[0], truth)])

        rows.append(dict(source=s, size=int(sizes[s]), n=len(tg),
                         prof=a.mean(), block=b.mean(), p=mcnemar(a, b),
                         gprof=ga.mean(), gblock=gb.mean(),
                         gp=mcnemar(ga, gb)))
        pooled["prof"].append(a); pooled["block"].append(b)
        pooled["gprof"].append(ga); pooled["gblock"].append(gb)
        pooled["obs"].append(tg.ksat_cmh.to_numpy(float))
        for i in range(3):
            pooled["kb"][i].append((p_prof[i + 1], p_block[i + 1]))

    print("EXACT USDA CLASS -- profile-level LOO vs whole-source hold-out")
    print(f"{'source':<22s} {'size':>6s} {'n':>5s} {'profile':>8s} "
          f"{'blocked':>8s} {'delta':>8s} {'p':>7s}")
    print("-" * 70)
    for r in sorted(rows, key=lambda r: -r["size"]):
        print(f"{r['source']:<22s} {r['size']:6d} {r['n']:5d} "
              f"{r['prof']*100:7.1f}% {r['block']*100:7.1f}% "
              f"{(r['block']-r['prof'])*100:+7.1f}pp {r['p']:7.3f}")
    A, B = np.concatenate(pooled["prof"]), np.concatenate(pooled["block"])
    print("-" * 70)
    print(f"{'POOLED':<22s} {len(gshp):6d} {len(A):5d} {A.mean()*100:7.1f}% "
          f"{B.mean()*100:7.1f}% {(B.mean()-A.mean())*100:+7.1f}pp "
          f"{mcnemar(A, B):7.3f}")

    print("\nAGGREGATED GROUP")
    print(f"{'source':<22s} {'size':>6s} {'n':>5s} {'profile':>8s} "
          f"{'blocked':>8s} {'delta':>8s} {'p':>7s}")
    print("-" * 70)
    for r in sorted(rows, key=lambda r: -r["size"]):
        print(f"{r['source']:<22s} {r['size']:6d} {r['n']:5d} "
              f"{r['gprof']*100:7.1f}% {r['gblock']*100:7.1f}% "
              f"{(r['gblock']-r['gprof'])*100:+7.1f}pp {r['gp']:7.3f}")
    GA, GB = np.concatenate(pooled["gprof"]), np.concatenate(pooled["gblock"])
    print("-" * 70)
    print(f"{'POOLED':<22s} {len(gshp):6d} {len(GA):5d} {GA.mean()*100:7.1f}% "
          f"{GB.mean()*100:7.1f}% {(GB.mean()-GA.mean())*100:+7.1f}pp "
          f"{mcnemar(GA, GB):7.3f}")

    obs = np.concatenate(pooled["obs"])
    km.report("KSAT", [
        ("profile LOO", km.ksat_scores(
            obs, *[np.concatenate([p[0] for p in pooled["kb"][i]])
                   for i in range(3)])),
        ("source-blocked", km.ksat_scores(
            obs, *[np.concatenate([p[1] for p in pooled["kb"][i]])
                   for i in range(3)])),
    ])


if __name__ == "__main__":
    main()
