"""Sensitivity of the vertisol texture result to the shrinkage cap.

Each soil is predicted as a new source (its own table removed from the
merged reference and from the classifier's training data), with the curve
converted from one density to a shrinking density BD(w) = min(cap,
1/(1/2.65 + w)), for:
  one density   the curve as one density gives it
  cap 1.4, 1.6  the values in the memo (1.6: Laikipia's published curves)
  cap 1.75      the median oven-dry clod density of KSSL Vertisols with clay
                45-65 % and density at 33 kPa below 1.25 (828 horizons)
  per sample    cap = 0.615 * BD1 + 1.011, the KSSL Vertisol regression of
                oven-dry on 33 kPa clod density (clay 45-65 %)
Laikipia: the published gravimetric curves (Willard et al. 2026, HydroShare,
doi:10.4211/hs.99080b48206e4b2c82e716516d463a61) with the KALRO density of the
core as BD1, texture from data/willard_reference.csv (70 samples, the
dispersed ones left out). The cap 1.6 arm is the published volumetric series,
the one the tool's reference holds. Tomejil: prepare_andalusia._tomejil_curves, points up to 1500 kPa,
one density = treatment mean. Run from the repository root.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
os.chdir(ROOT)
sys.argv = sys.argv[:1]

import swcc_texture as st                                    # noqa: E402
from verify_common import map_folds                          # noqa: E402

KPA_PER_CM = 0.0980665
ARMS = ["one density", "cap 1.4", "cap 1.6", "cap 1.75", "per sample"]


def bd_shrink(w, cap):
    return np.minimum(cap, 1.0 / (1.0 / 2.65 + np.asarray(w, float)))


def samples(soil):
    """[(id, h kPa, w, BD1, clay measured, class measured)]"""
    out = []
    if soil == "laikipia":
        # The published gravimetric curves (HydroShare, prepare_willard.RAW) at
        # the KALRO density of the core, the one the published volumetric
        # series starts from; its cap is 1.6.
        import prepare_willard as pw
        d = pw.read()
        p = d[d.SampleType == "point_sample"]
        ref = pd.read_csv(os.path.join("data", "willard_reference.csv")).set_index("layer_id")
        kalro = p[p.Measurement == "BulkDensity_KALRO"].set_index("SampleID").Value
        grav = p[p.Measurement.str.match(r"^pF[\d.]+_soilmoisture_gravimetric$")].dropna(subset=["Value"])
        for sid, x in grav.groupby("SampleID"):
            if f"WL_{sid}" not in ref.index:          # the dispersed samples
                continue
            pf = x.Measurement.str.extract(r"^pF([\d.]+)_")[0].astype(float).to_numpy()
            o = np.argsort(pf)
            h = np.where(pf[o] == 0, 0.0, 10.0 ** pf[o] * KPA_PER_CM)
            r = ref.loc[f"WL_{sid}"]
            out.append((sid, h, x.Value.to_numpy()[o] / 100.0, kalro[sid], r.clay, r.texture_class))
    else:
        sys.path.insert(0, ROOT)
        import prepare_andalusia as pa
        one = {"CT": 1.085, "DD": 1.125}
        for name, (t, _, h, w, _) in pa._tomejil_curves().items():
            keep = h * KPA_PER_CM <= 1500.0 * 1.04
            out.append((name, h[keep] * KPA_PER_CM, w[keep], one[t], pa.TOMEJIL_CLAY[t], "clay"))
    return out


def run(item, g):
    soil, arm = item
    src = "Willard_Laikipia" if soil == "laikipia" else "Vanderlinden_Tomejil"
    train = g[g.source_db != src].reset_index(drop=True)
    ref, clf = st.GshpReference(df=train), st.TextureGBM(df=train)
    res = []
    for sid, h, w, bd1, clay, cls in samples(soil):
        if arm == "one density":
            bd = np.full_like(w, bd1)
        elif arm == "per sample":
            bd = bd_shrink(w, 0.615 * bd1 + 1.011)
        else:
            bd = bd_shrink(w, float(arm.split()[1]))
        e = st.estimate(h, w * bd, ref=ref, clf=clf, sample_type="undisturbed")
        res.append(dict(soil=soil, arm=arm, sid=sid, cls=e["texture_class"], clay_p=e["fractions"]["clay"],
                        clay=clay, truth=cls, ks_knn=e["ksat"]["median_cmh"],
                        ks_phys=e["ksat"]["physical_matched_cmh"]))
    return res


def main():
    g = st.load_reference_df(st.DEFAULT_REFERENCE).reset_index(drop=True)
    g["source_db"] = g.source_db.fillna("KSSL")
    items = [(s, a) for s in ("laikipia", "tomejil") for a in ARMS]
    d = pd.DataFrame([r for part in map_folds(run, items, g, threads=2) for r in part])
    d.to_csv(os.path.join(HERE, "test_caps_results.csv"), index=False)
    for soil, x in d.groupby("soil", sort=False):
        print(f"\n{soil}: {x.sid.nunique()} samples, measured clay {x.groupby('sid').clay.first().mean():.1f} %")
        for arm in ARMS:
            y = x[x.arm == arm]
            cl = y[y.truth == "clay"]
            print(f"  {arm:12s} exact {np.mean(y.cls == y.truth) * 100:5.1f} %   clays read as clay "
                  f"{np.mean(cl.cls == 'clay') * 100:5.1f} % (n={len(cl)})   mean predicted clay "
                  f"{y.clay_p.mean():5.1f} %   clay error {np.mean(np.abs(y.clay_p - y.clay)):5.1f}   "
                  f"Ks kNN {y.ks_knn.median():.3g}, physical {y.ks_phys.median():.3g} cm/h, within 10x of each "
                  f"other {np.mean(np.abs(np.log10(y.ks_knn / y.ks_phys)) <= 1) * 100:.0f} %")
    # Laikipia's Ks is the Saturo infiltrometer (26 samples), the set's Ks in the reference
    lk = d[d.soil == "laikipia"].copy()
    ks = pd.read_csv(os.path.join("data", "willard_reference.csv")).set_index("layer_id").ksat_cmh
    lk["saturo"] = ("WL_" + lk.sid).map(ks)
    lk = lk.dropna(subset=["saturo"])
    print(f"\nlaikipia Ks against Saturo ({lk.sid.nunique()} samples, median {lk.groupby('sid').saturo.first().median():.1f} cm/h):")
    for arm in ARMS:
        y = lk[lk.arm == arm]
        for col, lab in (("ks_knn", "kNN"), ("ks_phys", "physical")):
            e = np.log10(y[col] / y.saturo)
            print(f"  {arm:12s} {lab:8s} median {y[col].median():.3g} cm/h  within 2x {np.mean(np.abs(e) <= np.log10(2)) * 100:3.0f} %"
                  f"  within 10x {np.mean(np.abs(e) <= 1) * 100:3.0f} %")


if __name__ == "__main__":
    main()
