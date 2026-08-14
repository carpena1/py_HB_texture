"""Independent verification of swcc_texture.py against the ROSETTA class-mean
van Genuchten parameters (Schaap, Leij & van Genuchten 2001, J. Hydrology
251:163-176; class-average lookup table = ROSETTA model H1).

For each of the 12 USDA classes we synthesize a noise-free retention curve
theta(psi) from the ROSETTA class-mean parameters (psi = 0..150 kPa, same
grid used for the Carsel benchmark), feed it to the tool, and check the
recovered class, particle fractions (vs USDA class centroid) and Ks.

ROSETTA table units: theta cm3/cm3; alpha as log10(1/cm); n as log10;
Ks as log10(cm/day). Converted here to the tool's conventions (alpha kPa^-1,
Ks cm/h) using 1 cm water head = 0.0980665 kPa and 24 h/day.
"""

import numpy as np

import swcc_texture as st

CM_TO_KPA = 0.0980665   # 1 cm water head in kPa

# class: (theta_r, theta_s, log10(alpha[1/cm]), log10(n), log10(Ks[cm/day]))
ROSETTA = {
    "sand":            (0.053, 0.375, -1.453, 0.502, 2.808),
    "loamy sand":      (0.049, 0.390, -1.459, 0.242, 2.022),
    "sandy loam":      (0.039, 0.387, -1.574, 0.161, 1.583),
    "loam":            (0.061, 0.399, -1.954, 0.168, 1.081),
    "silt":            (0.050, 0.489, -2.182, 0.225, 1.641),
    "silt loam":       (0.065, 0.439, -2.296, 0.221, 1.261),
    "sandy clay loam": (0.063, 0.384, -1.676, 0.124, 1.120),
    "clay loam":       (0.079, 0.442, -1.801, 0.151, 0.913),
    "silty clay loam": (0.090, 0.482, -2.076, 0.182, 1.046),
    "sandy clay":      (0.117, 0.385, -1.476, 0.082, 1.055),
    "silty clay":      (0.111, 0.481, -1.790, 0.121, 0.983),
    "clay":            (0.098, 0.459, -1.825, 0.098, 1.169),
}
ORDER = list(ROSETTA)


def rosetta_params(cls):
    tr, ts, la, ln, lks = ROSETTA[cls]
    alpha_kpa = (10.0 ** la) / CM_TO_KPA
    n = 10.0 ** ln
    ks_cmh = (10.0 ** lks) / 24.0
    return tr, ts, alpha_kpa, n, ks_cmh


def main():
    centroids = st.usda_centroids()
    ref = st.GshpReference()
    clf = st.TextureGBM(df=st.load_reference_df())
    h = np.arange(0.0, 150.0 + 0.1, 5.0)

    n_top1 = n_top2 = n_ks_in = n_frac_in = 0
    frac_err, ks_logerr = [], []

    print(f"{'true class':<16s} {'predicted':<16s} {'P(true)':>7s} "
          f"{'rank':>4s}  {'sand/silt/clay pred':<20s} {'centroid':<15s} "
          f"{'Ks pred':>8s} {'Ks true':>8s} {'nb':>4s} {'in5-95%':>7s}")
    for name in ORDER:
        tr, ts, alpha, n, ks_true = rosetta_params(name)
        theta = st.vg_theta(h, tr, ts, alpha, n)
        res = st.estimate(h, theta, ref=ref, clf=clf)

        probs = res["class_probabilities"]
        pred = res["texture_class"]
        p_true = probs.get(name, 0.0)
        rank = list(probs).index(name) + 1 if name in probs else len(ORDER)
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
        in_ks = ks["p5_cmh"] <= ks_true <= ks["p95_cmh"]
        n_ks_in += in_ks
        ks_logerr.append(abs(np.log10(ks["median_cmh"] / ks_true)))

        print(f"{name:<16s} {pred:<16s} {p_true:7.2f} {rank:4d}  "
              f"{fr['sand']:4.0f}/{fr['silt']:4.0f}/{fr['clay']:4.0f}{'':>7s}"
              f"{cen[0]:3.0f}/{cen[1]:3.0f}/{cen[2]:3.0f}{'':>4s}"
              f"{ks['median_cmh']:8.3g} {ks_true:8.3g} "
              f"{ks['n_neighbors_with_ksat']:4.0f} "
              f"{'yes' if in_ks else 'NO':>7s}")

    n = len(ORDER)
    print(f"\nSummary over {n} classes (ROSETTA class-mean benchmark):")
    print(f"  exact class match:            {n_top1}/{n}")
    print(f"  true class in top 2:          {n_top2}/{n}")
    print(f"  mean |fraction error| vs centroid: {np.mean(frac_err):.1f} % "
          f"(centroid inside 5-95 % range: {n_frac_in}/{n})")
    print(f"  Ks: mean |log10 error| = {np.mean(ks_logerr):.2f} "
          f"(typical factor {10 ** np.mean(ks_logerr):.1f}), "
          f"true Ks inside 5-95 % range: {n_ks_in}/{n}")


if __name__ == "__main__":
    main()
