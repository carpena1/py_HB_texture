"""Leave-one-out verification of swcc_texture.py against the source GSHP
database, using 12 representative soils (one per USDA class).

For each class we take the GSHP layer closest to that class's median in the
tool's van Genuchten feature space (preferring layers that carry a measured
Ksat) -- the same soils shipped as Testing/<i>_<code>_GHSP.csv. We synthesize
its retention curve, then classify it with that layer REMOVED from the
reference database (leave-one-out), so the soil cannot match itself.

Unlike the Carsel/ROSETTA benchmarks, these are real soils with measured
texture and Ksat, so predictions are checked against the actual per-soil
sand/silt/clay and Ks -- not class centroids.
"""

import numpy as np
import pandas as pd

import swcc_texture as st

ORDER = ["sand", "loamy sand", "sandy loam", "loam", "silt", "silt loam",
         "sandy clay loam", "clay loam", "silty clay loam", "sandy clay",
         "silty clay", "clay"]


def feat(d):
    return np.column_stack([np.log10(d.alpha_kpa), np.log10(d.n - 1.0),
                            d.thetar, d.thetas])


def representative_layer(df, cls):
    """Index of the layer closest to the class median in feature space,
    preferring layers with a measured Ksat (matches the sample-file generator)."""
    g = df[df.texture_class == cls]
    pool = g[g.ksat_cmh.notna()]
    if len(pool) == 0:
        pool = g
    med = np.median(feat(g), axis=0)
    z = (feat(pool) - med) / (feat(g).std(axis=0) + 1e-9)
    return pool.index[np.argmin((z ** 2).sum(1))]


def main():
    df = pd.read_csv(st.REFERENCE_CSV)
    h = np.arange(0.0, 150.0 + 0.1, 5.0)

    n_top1 = n_top2 = n_ks_in = n_frac_in = 0
    frac_err, ks_logerr = [], []
    ks_n = 0

    print(f"{'true class':<16s} {'predicted':<16s} {'P(true)':>7s} {'rank':>4s}  "
          f"{'s/si/cl pred':<14s} {'s/si/cl true':<14s} "
          f"{'Ks pred':>8s} {'Ks true':>8s} {'nb':>4s} {'Ks in5-95%':>10s}")
    for name in ORDER:
        j = representative_layer(df, name)
        row = df.loc[j]
        theta = st.vg_theta(h, row.thetar, row.thetas, row.alpha_kpa, row.n)

        # leave-one-out: drop this exact layer from the reference
        loo = df.drop(index=j).reset_index(drop=True)
        ref = st.GshpReference(df=loo)
        res = st.estimate(h, theta, ref=ref)

        probs = res["class_probabilities"]
        pred = res["texture_class"]
        p_true = probs.get(name, 0.0)
        rank = list(probs).index(name) + 1 if name in probs else len(ORDER)
        n_top1 += pred == name
        n_top2 += rank <= 2

        fr = res["fractions"]
        true_frac = (row["sand"], row["silt"], row["clay"])
        frac_err.append(np.mean([abs(fr[c] - tv) for c, tv
                                 in zip(("sand", "silt", "clay"), true_frac)]))
        in_frac = all(fr["p5"][c] <= tv <= fr["p95"][c]
                      for c, tv in zip(("sand", "silt", "clay"), true_frac))
        n_frac_in += in_frac

        ks = res["ksat"]
        has_ks = np.isfinite(row.ksat_cmh)
        in_ks = has_ks and ks["p5_cmh"] <= row.ksat_cmh <= ks["p95_cmh"]
        if has_ks:
            ks_n += 1
            n_ks_in += in_ks
            ks_logerr.append(abs(np.log10(ks["median_cmh"] / row.ksat_cmh)))
        ks_true_s = f"{row.ksat_cmh:8.3g}" if has_ks else f"{'NA':>8s}"

        print(f"{name:<16s} {pred:<16s} {p_true:7.2f} {rank:4d}  "
              f"{fr['sand']:3.0f}/{fr['silt']:3.0f}/{fr['clay']:3.0f}{'':>4s}"
              f"{true_frac[0]:3.0f}/{true_frac[1]:3.0f}/{true_frac[2]:3.0f}{'':>4s}"
              f"{ks['median_cmh']:8.3g} {ks_true_s} "
              f"{ks['n_neighbors_with_ksat']:4.0f} "
              f"{('yes' if in_ks else 'NO') if has_ks else '-':>10s}")

    n = len(ORDER)
    print(f"\nLeave-one-out over 12 representative GSHP soils (target soil "
          f"removed from the reference before each prediction):")
    print(f"  exact class match:            {n_top1}/{n}")
    print(f"  true class in top 2:          {n_top2}/{n}")
    print(f"  mean |fraction error| vs measured: {np.mean(frac_err):.1f} % "
          f"(measured inside 5-95 % range: {n_frac_in}/{n})")
    print(f"  Ks (of {ks_n} soils with measured Ksat): "
          f"mean |log10 error| = {np.mean(ks_logerr):.2f} "
          f"(typical factor {10 ** np.mean(ks_logerr):.1f}), "
          f"measured Ks inside 5-95 % range: {n_ks_in}/{ks_n}")


if __name__ == "__main__":
    main()
