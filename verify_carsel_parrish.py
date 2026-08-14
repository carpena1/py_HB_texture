"""Independent verification of swcc_texture.py against the Carsel & Parrish
(1988) class-typical benchmark.

The benchmark (used ONLY here, never during development) holds, for each of
the 12 USDA texture classes, a synthetic retention curve theta(psi) at
psi = 0..150 kPa. Curves are read from the wide file
Testing_Carsel_Parrish/Carsel_Parrish.csv; the true Ks (cm/h) values are the
Carsel & Parrish class-typical constants.

Checks per class:
  - predicted USDA class vs true class (and probability given to truth)
  - predicted mean sand/silt/clay vs the centroid of the true class in the
    USDA texture triangle (per agreed protocol), incl. 5-95 % range coverage
  - predicted Ks vs benchmark Ks, incl. 5-95 % range coverage
"""

import csv

import numpy as np

import swcc_texture as st

BENCHMARK_CSV = "data/Carsel_Parrish.csv"

# short code (row 2 of the wide file) -> canonical USDA class name
CODE_TO_NAME = {
    "S": "sand", "LS": "loamy sand", "SL": "sandy loam", "L": "loam",
    "Si": "silt", "SiL": "silt loam", "SCL": "sandy clay loam",
    "CL": "clay loam", "SiCL": "silty clay loam", "SC": "sandy clay",
    "SiC": "silty clay", "C": "clay",
}

# Carsel & Parrish (1988) class-typical Ks in cm/h.
KS_TRUE = {
    "sand": 29.7, "loamy sand": 14.5917, "sandy loam": 4.420833,
    "loam": 1.04, "silt": 0.25, "silt loam": 0.45,
    "sandy clay loam": 1.31, "clay loam": 0.26, "silty clay loam": 0.07,
    "sandy clay": 0.12, "silty clay": 0.02, "clay": 0.2,
}


def load_benchmark():
    rows = list(csv.reader(open(BENCHMARK_CSV, encoding="utf-8-sig")))
    codes = rows[1]
    names = [CODE_TO_NAME[codes[2 * (i - 1) + 1].strip()] for i in range(1, 13)]

    curves = {name: {"h": [], "theta": []} for name in names}
    for r in rows[3:]:
        for i, name in enumerate(names, start=1):
            sc = 2 * (i - 1)
            if r[sc].strip() == "":
                continue
            curves[name]["h"].append(float(r[sc]))
            curves[name]["theta"].append(float(r[sc + 1]))
    return names, curves, KS_TRUE


def reference_counts():
    """Per-class curated GSHP counts: n_ref (texture/fraction layers) and
    n_ks (subset with a measured Ksat)."""
    import pandas as pd
    df = st.load_reference_df()
    n_ref = df["texture_class"].value_counts()
    n_ks = df[df["ksat_cmh"].notna()]["texture_class"].value_counts()
    return n_ref, n_ks


def main():
    names, curves, ks_true = load_benchmark()
    centroids = st.usda_centroids()
    ref = st.GshpReference()
    clf = st.TextureGBM(df=st.load_reference_df())
    n_ref, n_ks = reference_counts()

    n_top1 = n_top2 = n_ks_in = n_frac_in = 0
    frac_err, ks_logerr = [], []

    print(f"{'true class':<16s} {'n_ref':>6s} {'n_ks':>5s} | "
          f"{'predicted':<16s} {'P(true)':>7s} {'rank':>4s}  "
          f"{'sand/silt/clay pred':<20s} {'centroid':<15s} "
          f"{'Ks pred':>8s} {'Ks true':>8s} {'nb':>4s} {'in5-95%':>7s}")
    for name in names:
        cur = curves[name]
        res = st.estimate(cur["h"], cur["theta"], ref=ref, clf=clf)
        probs = res["class_probabilities"]
        pred = res["texture_class"]
        p_true = probs.get(name, 0.0)
        rank = list(probs).index(name) + 1 if name in probs else len(st.USDA_CLASSES)
        n_top1 += pred == name
        n_top2 += rank <= 2

        fr = res["fractions"]
        cen = centroids[name]
        frac_err.append(np.mean([abs(fr[c] - cv) for c, cv
                                 in zip(("sand", "silt", "clay"), cen)]))
        in_frac = all(fr["p5"][c] <= cv <= fr["p95"][c]
                      for c, cv in zip(("sand", "silt", "clay"), cen))
        n_frac_in += in_frac

        ks = res["ksat"]
        in_ks = ks["p5_cmh"] <= ks_true[name] <= ks["p95_cmh"]
        n_ks_in += in_ks
        ks_logerr.append(abs(np.log10(ks["median_cmh"] / ks_true[name])))

        print(f"{name:<16s} {n_ref[name]:6d} {n_ks[name]:5d} | "
              f"{pred:<16s} {p_true:7.2f} {rank:4d}  "
              f"{fr['sand']:4.0f}/{fr['silt']:4.0f}/{fr['clay']:4.0f}{'':>7s}"
              f"{cen[0]:3.0f}/{cen[1]:3.0f}/{cen[2]:3.0f}{'':>4s}"
              f"{ks['median_cmh']:8.3g} {ks_true[name]:8.3g} "
              f"{ks['n_neighbors_with_ksat']:4.0f} "
              f"{'yes' if in_ks else 'NO':>7s}")

    n = len(names)
    print(f"\nReference database (curated): {int(n_ref.sum())} layers for "
          f"texture/fractions, {int(n_ks.sum())} with measured Ksat.")
    print(f"('nb' = mean number of the k=30 nearest neighbors that carried a "
          f"measured Ksat, i.e. the support behind each Ks estimate.)")
    print(f"\nSummary over {n} classes:")
    print(f"  exact class match:            {n_top1}/{n}")
    print(f"  true class in top 2:          {n_top2}/{n}")
    print(f"  mean |fraction error| vs centroid: {np.mean(frac_err):.1f} % "
          f"(centroid inside 5-95 % range: {n_frac_in}/{n})")
    print(f"  Ks: mean |log10 error| = {np.mean(ks_logerr):.2f} "
          f"(i.e. typical factor {10 ** np.mean(ks_logerr):.1f}), "
          f"true Ks inside 5-95 % range: {n_ks_in}/{n}")


if __name__ == "__main__":
    main()
