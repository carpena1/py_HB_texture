"""Accuracy at the level of aggregated texture groups.

The 12 USDA classes are collapsed into four coarse groups; each verification's
per-class prediction is mapped to its group, so a miss that stays inside the
same group counts as correct. Reuses the three verification drivers.

Four overarching groups (Sandy / Loamy / Silty / Clayey):
  Sandy : sand, loamy sand
  Loamy : sandy loam, loam, sandy clay loam, clay loam
  Silty : silt, silt loam, silty clay loam
  Clayey: sandy clay, silty clay, clay
"""

import numpy as np
import pandas as pd

import swcc_texture as st
import verify_carsel_parrish as vcp
import verify_rosetta as vr
import verify_gshp as vg

GROUP = {
    "sand": "Sandy", "loamy sand": "Sandy",
    "sandy loam": "Loamy", "loam": "Loamy",
    "sandy clay loam": "Loamy", "clay loam": "Loamy",
    "silt": "Silty", "silt loam": "Silty", "silty clay loam": "Silty",
    "sandy clay": "Clayey", "silty clay": "Clayey", "clay": "Clayey",
}
GROUP_ORDER = ["Sandy", "Loamy", "Silty", "Clayey"]

H = np.arange(0.0, 150.0 + 0.1, 5.0)


def preds_carsel():
    names, curves, _ = vcp.load_benchmark()
    ref = st.GshpReference()
    return [(n, st.estimate(curves[n]["h"], curves[n]["theta"],
                            ref=ref)["texture_class"]) for n in names]


def preds_rosetta():
    ref = st.GshpReference()
    out = []
    for n in vr.ORDER:
        tr, ts, a, nn, _ = vr.rosetta_params(n)
        theta = st.vg_theta(H, tr, ts, a, nn)
        out.append((n, st.estimate(H, theta, ref=ref)["texture_class"]))
    return out


def preds_gshp():
    df = pd.read_csv(st.REFERENCE_CSV)
    out = []
    for n in vg.ORDER:
        j = vg.representative_layer(df, n)
        row = df.loc[j]
        theta = st.vg_theta(H, row.thetar, row.thetas, row.alpha_kpa, row.n)
        loo = st.GshpReference(df=df.drop(index=j).reset_index(drop=True))
        out.append((n, st.estimate(H, theta, ref=loo)["texture_class"]))
    return out


def report(title, preds):
    n = len(preds)
    n_class = sum(t == p for t, p in preds)
    n_group = sum(GROUP[t] == GROUP[p] for t, p in preds)
    print(f"\n=== {title} ===")
    print(f"  class-level:  {n_class}/{n} exact")
    print(f"  group-level:  {n_group}/{n} correct group")
    misses = [(t, p) for t, p in preds if GROUP[t] != GROUP[p]]
    if misses:
        print("  cross-group misses (true -> pred | group true -> group pred):")
        for t, p in misses:
            print(f"    {t:<16s} -> {p:<16s} | {GROUP[t]:<24s} -> {GROUP[p]}")
    return n_group, n


def main():
    print("Aggregated texture-group accuracy\n" + "=" * 34)
    print("Groups:")
    for g in GROUP_ORDER:
        members = [c for c, gg in GROUP.items() if gg == g]
        print(f"  {g:<24s}: {', '.join(members)}")

    totals = []
    for title, fn in [("Carsel & Parrish (1988)", preds_carsel),
                      ("ROSETTA class means", preds_rosetta),
                      ("GSHP leave-one-out", preds_gshp)]:
        ng, n = report(title, fn())
        totals.append((title, ng, n))

    print("\n=== Summary (group-level) ===")
    for title, ng, n in totals:
        print(f"  {title:<26s} {ng}/{n}  ({100*ng/n:.0f} %)")


if __name__ == "__main__":
    main()
