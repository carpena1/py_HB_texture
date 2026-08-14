"""Verification against the HYPRES class pedotransfer functions
(Wosten et al. 1999; ESDAC / European Soil Database v2.0 distribution).

WHY THIS IS A BENCHMARK AND NOT A REFERENCE
    The plan was to merge HYPRES's ~4,486 European horizons into the kNN
    reference. That is not possible: the released HYPRES v1.0 contains only
    metadata plus the class and continuous pedotransfer functions. Its
    HYPRES_Readme states the release "does not contain ... HYPRES project
    source data and results as no agreement has been reached with the
    participating institutions regarding their distribution."  So the
    sample-level tables (RAWRET / RAWK / HYDRAULIC_PROPS) were never released.

    What we do have is 11 class-average Mualem-van Genuchten parameter sets.
    That is far too few to matter as reference points among ~10,000, but it
    makes a genuine *European* benchmark -- a counterpart to the US-centric
    Carsel & Parrish and ROSETTA benchmarks.

HYPRES uses FAO/SGDBE texture classes (coarse .. very fine, split into topsoil
and subsoil), not the 12 USDA classes, so the tool's USDA prediction is scored
by whether its predicted sand/silt/clay falls inside the correct FAO class
envelope. Those envelopes (percent of the < 2 mm fraction) are:

    Coarse       clay < 18 and sand > 65
    Medium       (18 <= clay < 35 and sand >= 15) or
                 (clay < 18 and 15 <= sand <= 65)
    Medium fine  clay < 35 and sand < 15
    Fine         35 <= clay < 60
    Very fine    clay >= 60

The organic class is excluded: it is not a mineral texture class.

Parameter units as distributed: alpha in 1/cm, Ks in cm/day.
"""

import numpy as np

import swcc_texture as st

CM_TO_KPA = 0.0980665

# class -> (theta_r, theta_s, alpha[1/cm], n, m, Ks[cm/day])
HYPRES = {
    ("topsoil", "coarse"):      (0.025, 0.403, 0.0383, 1.3774, 0.2740, 60.000),
    ("topsoil", "medium"):      (0.010, 0.439, 0.0314, 1.1804, 0.1528, 12.061),
    ("topsoil", "medium fine"): (0.010, 0.430, 0.0083, 1.2539, 0.2025,  2.272),
    ("topsoil", "fine"):        (0.010, 0.520, 0.0367, 1.1012, 0.0919, 24.800),
    ("topsoil", "very fine"):   (0.010, 0.614, 0.0265, 1.1033, 0.0936, 15.000),
    ("subsoil", "coarse"):      (0.025, 0.366, 0.0430, 1.5206, 0.3424, 70.000),
    ("subsoil", "medium"):      (0.010, 0.392, 0.0249, 1.1689, 0.1445, 10.755),
    ("subsoil", "medium fine"): (0.010, 0.412, 0.0082, 1.2179, 0.1789,  4.000),
    ("subsoil", "fine"):        (0.010, 0.481, 0.0198, 1.0861, 0.0793,  8.500),
    ("subsoil", "very fine"):   (0.010, 0.538, 0.0168, 1.0730, 0.0680,  8.235),
}

# Representative depth for each pedological position, used with --depth.
DEPTH_CM = {"topsoil": 15.0, "subsoil": 60.0}


def fao_class(sand, silt, clay):
    if clay < 18 and sand > 65:
        return "coarse"
    if (18 <= clay < 35 and sand >= 15) or (clay < 18 and 15 <= sand <= 65):
        return "medium"
    if clay < 35 and sand < 15:
        return "medium fine"
    if 35 <= clay < 60:
        return "fine"
    return "very fine"


def main():
    h = np.arange(0.0, 150.0 + 0.1, 5.0)
    ref_plain = st.GshpReference(reference="gshp")

    import pandas as pd
    ref_depth = st.GshpReference(df=pd.read_csv(st.REFERENCE_CSV),
                                 use_depth=True)

    n_ok = n_ok_d = 0
    ks_log, ks_in = [], 0
    print(f"{'HYPRES class':<22s} {'pred USDA':<16s} "
          f"{'s/si/cl pred':<14s} {'FAO pred':<12s} {'ok':>3s} "
          f"{'+depth':>7s} {'Ks pred':>8s} {'Ks HYPRES':>10s}")
    for (pos, cls), (tr, ts, a_cm, n, m, ks_day) in HYPRES.items():
        # The distribution reports m, which should equal 1 - 1/n (Mualem).
        assert abs(m - (1 - 1 / n)) < 5e-3, (cls, m, 1 - 1 / n)
        alpha = a_cm / CM_TO_KPA
        ks_true = ks_day / 24.0
        theta = st.vg_theta(h, tr, ts, alpha, n)

        res = st.estimate(h, theta, ref=ref_plain)
        fr = res["fractions"]
        pred_fao = fao_class(fr["sand"], fr["silt"], fr["clay"])
        ok = pred_fao == cls
        n_ok += ok

        rd = st.estimate(h, theta, ref=ref_depth, depth=DEPTH_CM[pos])
        fd = rd["fractions"]
        ok_d = fao_class(fd["sand"], fd["silt"], fd["clay"]) == cls
        n_ok_d += ok_d

        ks = res["ksat"]
        if np.isfinite(ks["median_cmh"]):
            ks_log.append(abs(np.log10(ks["median_cmh"] / ks_true)))
            ks_in += ks["p5_cmh"] <= ks_true <= ks["p95_cmh"]

        print(f"{pos+' '+cls:<22s} {res['texture_class']:<16s} "
              f"{fr['sand']:4.0f}/{fr['silt']:3.0f}/{fr['clay']:3.0f}{'':>2s}"
              f"{pred_fao:<12s} {'yes' if ok else 'NO':>3s} "
              f"{'yes' if ok_d else 'NO':>7s} "
              f"{ks['median_cmh']:8.3g} {ks_true:10.3g}")

    n = len(HYPRES)
    print(f"\nFAO texture class correct: {n_ok}/{n}   "
          f"with --depth: {n_ok_d}/{n}")
    print(f"Ks: mean |log10 error| {np.mean(ks_log):.2f} "
          f"(typical factor {10**np.mean(ks_log):.1f}), "
          f"HYPRES Ks inside 5-95 % range {ks_in}/{len(ks_log)}")


if __name__ == "__main__":
    main()
