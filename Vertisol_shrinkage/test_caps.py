"""Sensitivity of the vertisol texture result to the shrinkage cap.

Each soil is predicted as a new source (its own table removed from the
merged reference and from the classifier's training data), with the curve
converted from one density to a shrinking density BD(w) = min(cap,
1/(1/2.65 + w)), for:
  one density   the curve as one density gives it
  cap 1.4, 1.6  the values in the memo
  cap 1.75      the median oven-dry clod density of KSSL Vertisols with clay
                45-65 % and density at 33 kPa below 1.25 (828 horizons)
  per sample    cap = 0.615 * BD1 + 1.011, the KSSL Vertisol regression of
                oven-dry on 33 kPa clod density (clay 45-65 %)
Laikipia: the one-density curve of prepare_willard.curve (pF 0 at 95 % of
porosity, BD1 from porosity; prepare_willard.scaled_curve) before its shrinkage step, with texture from
data/willard_reference.csv (86 samples), so the cap 1.75 arm is the curve the
tool's reference holds. Tomejil: prepare_andalusia._tomejil_curves, points up to 1500 kPa,
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
        import prepare_willard as pw
        ref = pd.read_csv(os.path.join("data", "willard_reference.csv")).set_index("layer_id")
        cap = pw.SHRINK_CAP
        pw.SHRINK_CAP = np.inf                       # no cap: th = w / (1/2.65 + w)
        try:
            curves = pw.scaled_curve(pw.read())
        finally:
            pw.SHRINK_CAP = cap
        por = pw.porosity(pw.read())
        for sid, (h, th, _) in curves.items():
            if f"WL_{sid}" not in ref.index:          # left out of the reference
                continue
            w = th / (pw.PARTICLE_DENSITY * (1.0 - th))   # the gravimetric curve back
            r = ref.loc[f"WL_{sid}"]
            out.append((sid, h, w, pw.PARTICLE_DENSITY * (1.0 - por[sid]), r.clay, r.texture_class))
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


if __name__ == "__main__":
    main()
