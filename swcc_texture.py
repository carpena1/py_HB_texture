"""Estimate USDA soil texture and Ksat from a measured soil water
characteristic curve (SWCC).

Framework (see Problem_description.docx):
 1. Fit van Genuchten (vG) parameters (thetar, thetas, alpha, n) to the
    measured (h, theta) data, with the Mualem constraint m = 1 - 1/n.
 2. Infer USDA texture class, particle fractions (% sand/silt/clay) and
    saturated hydraulic conductivity Ks (cm/h) by distance-weighted
    k-nearest-neighbor lookup in the GSHP database (Gupta et al. 2022),
    using the fitted vG parameters as coordinates.

Uncertainty: vG fit uncertainty is propagated by Monte Carlo sampling of
the fit covariance; each draw votes through its own kNN neighborhood.

Input units: h in kPa (suction, positive); theta as volumetric fraction
(m3/m3) or percent — auto-detected (any theta > 1 means percent).

CLI:  swcc_texture.py datafile
      datafile: text/CSV file with two columns: h_kPa, theta
"""

import argparse
import json
import os
import sys

import numpy as np
from scipy.optimize import curve_fit

# Resolve reference tables relative to this script, so the tool works from any
# working directory.
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

# GSHP alone. Kept as a module constant because the historical verification
# scripts (Carsel & Parrish, ROSETTA, region, variants) were all scored against
# GSHP and their published tables must stay reproducible.
REFERENCE_CSV = os.path.join(DATA_DIR, "gshp_reference.csv")

# Named reference sets. "merged" is the default: GSHP plus the NCSS/KSSL
# layers, whose class mix complements GSHP's (1.8x the silty clay loam, only
# 4 % as much sand). The merge is not a significant accuracy win on its own
# (+0.7 pp, p=0.35) but it is not harmful either, it lifts the Loamy group
# past the 50 % bar under Bonferroni (p=0.019 -> 0.004), and it improves Ksat
# within a factor of two from 44 % to 47 % by changing which GSHP neighbours
# are selected -- despite KSSL carrying no Ksat of its own.
REFERENCE_SETS = {
    "merged": ["gshp_reference.csv", "kssl_reference.csv"],
    "gshp": ["gshp_reference.csv"],
    "kssl": ["kssl_reference.csv"],
    "all": ["gshp_reference.csv", "unsoda_reference.csv",
            "kssl_reference.csv"],
}
DEFAULT_REFERENCE = "merged"


def load_reference_df(name=DEFAULT_REFERENCE):
    """Concatenate the CSVs making up a named reference set."""
    import pandas as pd
    if name not in REFERENCE_SETS:
        raise ValueError(f"unknown reference set {name!r}; "
                         f"choose from {sorted(REFERENCE_SETS)}")
    return pd.concat([pd.read_csv(os.path.join(DATA_DIR, f))
                      for f in REFERENCE_SETS[name]], ignore_index=True)

USDA_CLASSES = [
    "sand", "loamy sand", "sandy loam", "loam", "silt", "silt loam",
    "sandy clay loam", "clay loam", "silty clay loam", "sandy clay",
    "silty clay", "clay",
]

# ---------------------------------------------------------------------------
# van Genuchten model and fitting
# ---------------------------------------------------------------------------

def vg_theta(h, thetar, thetas, alpha, n):
    """van Genuchten retention, m = 1 - 1/n, h in kPa >= 0."""
    m = 1.0 - 1.0 / n
    return thetar + (thetas - thetar) * (1.0 + (alpha * np.asarray(h)) ** n) ** (-m)


# Matric potentials (kPa) at which retention is sampled for curve-space
# matching: log-spaced from near saturation to the wilting point, and
# including the two classical anchors, field capacity (33) and 1500.
CURVE_HEADS = np.array([1.0, 3.0, 10.0, 33.0, 100.0, 330.0, 1000.0, 1500.0])


def _curve_features(thetar, thetas, alpha, n):
    """Retention sampled at CURVE_HEADS, shape (n_rows, n_heads).

    Matching on the curve rather than on (alpha, n) sidesteps the fit ridge:
    alpha and n trade off against each other, so two soils with nearly
    identical curves can sit far apart in parameter space, and two soils that
    are close in parameter space can have visibly different curves.
    """
    thetar, thetas, alpha, n = (np.atleast_1d(np.asarray(x, dtype=float))
                                for x in (thetar, thetas, alpha, n))
    m = 1.0 - 1.0 / n
    ah = alpha[:, None] * CURVE_HEADS[None, :]
    return (thetar[:, None] + (thetas - thetar)[:, None]
            * (1.0 + ah ** n[:, None]) ** (-m[:, None]))


def _vg_transformed(h, thetar, thetas, la, ln1):
    return vg_theta(h, thetar, thetas, 10.0 ** la, 1.0 + 10.0 ** ln1)


def fit_vg(h, theta):
    """Fit vG parameters to (h [kPa], theta [m3/m3]) data.

    Fits in transformed space p = [thetar, thetas, log10(alpha), log10(n-1)]
    for stability. Returns (popt, pcov) in that transformed space.
    """
    h = np.asarray(h, dtype=float)
    theta = np.asarray(theta, dtype=float)
    if len(h) < 5:
        raise ValueError("need at least 5 (h, theta) points to fit 4 parameters")

    tmin, tmax = theta.min(), theta.max()
    hpos = h[h > 0]
    p0 = [max(0.5 * tmin, 1e-3), tmax, np.log10(1.0 / np.median(hpos)), np.log10(0.5)]
    lb = [0.0, 0.5 * tmax, -4.0, np.log10(0.01)]
    ub = [tmin + 1e-9, 1.0, 1.5, np.log10(10.0)]
    p0 = np.clip(p0, lb, ub)

    popt, pcov = curve_fit(_vg_transformed, h, theta, p0=p0,
                           bounds=(lb, ub), maxfev=20000)
    return popt, pcov


# ---------------------------------------------------------------------------
# USDA texture triangle
# ---------------------------------------------------------------------------

def usda_class(sand, silt, clay):
    """USDA texture class from fractions in % (must sum to ~100)."""
    if silt + 1.5 * clay < 15:
        return "sand"
    if silt + 2.0 * clay < 30:
        return "loamy sand"
    if (7 <= clay <= 20 and sand > 52 and silt + 2 * clay >= 30) or \
       (clay < 7 and silt < 50 and silt + 2 * clay >= 30):
        return "sandy loam"
    if 7 <= clay <= 27 and 28 <= silt < 50 and sand <= 52:
        return "loam"
    if silt >= 80 and clay < 12:
        return "silt"
    if (silt >= 50 and 12 <= clay < 27) or (50 <= silt < 80 and clay < 12):
        return "silt loam"
    if 20 <= clay < 35 and silt < 28 and sand > 45:
        return "sandy clay loam"
    if 27 <= clay < 40 and 20 < sand <= 45:
        return "clay loam"
    if 27 <= clay < 40 and sand <= 20:
        return "silty clay loam"
    if clay >= 35 and sand > 45:
        return "sandy clay"
    if clay >= 40 and silt >= 40:
        return "silty clay"
    if clay >= 40 and sand <= 45 and silt < 40:
        return "clay"
    return "loam"  # boundary fall-through (rounding artifacts)


def usda_centroids(step=0.25):
    """Centroid (mean sand/silt/clay, %) of each USDA class polygon,
    computed on a regular grid over the texture triangle."""
    grid = np.arange(0.0, 100.0 + step / 2, step)
    sums = {c: np.zeros(4) for c in USDA_CLASSES}
    for sa in grid:
        for cl in np.arange(0.0, 100.0 - sa + step / 2, step):
            si = 100.0 - sa - cl
            c = usda_class(sa, si, cl)
            sums[c] += (sa, si, cl, 1.0)
    return {c: tuple(v[:3] / v[3]) for c, v in sums.items()}


# ---------------------------------------------------------------------------
# GSHP kNN inference
# ---------------------------------------------------------------------------

class GshpReference:
    def __init__(self, path=None, df=None, reference=DEFAULT_REFERENCE,
                 use_depth=False, tau=1.0, feature_mode="vg",
                 quality_weight=False):
        """Build the reference from the CSV at `path`, or from a preloaded
        DataFrame `df` (used for leave-one-out verification, where the target
        soil is dropped before constructing the reference).

        `use_depth=True` adds the sample mid-depth as a fifth, equally
        weighted feature (log10(1+depth_cm), standardized); reference rows
        without a depth are dropped in that case.
        """
        import pandas as pd
        if df is None:
            df = (pd.read_csv(path) if path is not None
                  else load_reference_df(reference))
        self.use_depth = use_depth
        if use_depth:
            df = df[df["depth_cm"].notna()].reset_index(drop=True)
        self.layer_id = df["layer_id"].to_numpy()
        self.profile_id = (df["profile_id"].to_numpy()
                           if "profile_id" in df.columns else None)
        self.classes = df["texture_class"].to_numpy()
        # Class-frequency vote weighting, w = (1/count)^tau.
        #
        # tau = 1 imposes a uniform prior over the 12 USDA classes instead of
        # GSHP's own distribution (sand alone is ~39 % of the database). That
        # maximises macro-recall but is brutal for rare classes: silt has 32
        # layers against sand's 3,930, so a single silt neighbour outvotes 123
        # sand neighbours and silt gets predicted 2.7x more often than it
        # occurs (precision 18.6 %). tau = 0 disables the correction entirely
        # and lets sand dominate. Intermediate values trade recall against
        # precision; see verify_tau.py.
        self.tau = tau
        counts = df["texture_class"].value_counts()
        self.class_weight = ((1.0 / counts) ** tau)[df["texture_class"]].to_numpy()
        self.feature_mode = feature_mode

        # Optional per-row reliability weight. GSHP publishes the standard
        # error of its own vG fit; rse_alpha/rse_n are those as relative
        # errors. A layer whose alpha is barely identified should not vote as
        # loudly as one measured over a full curve. Rows lacking the columns
        # (other databases) or the values keep weight 1.
        self.row_weight = np.ones(len(df))
        if quality_weight:
            ra = (pd.to_numeric(df["rse_alpha"], errors="coerce").to_numpy()
                  if "rse_alpha" in df.columns else np.zeros(len(df)))
            rn = (pd.to_numeric(df["rse_n"], errors="coerce").to_numpy()
                  if "rse_n" in df.columns else np.zeros(len(df)))
            ra = np.nan_to_num(ra, nan=0.0)
            rn = np.nan_to_num(rn, nan=0.0)
            self.row_weight = 1.0 / (1.0 + ra ** 2 + rn ** 2)
        self.fractions = df[["sand", "silt", "clay"]].to_numpy()
        self.ksat = df["ksat_cmh"].to_numpy()
        if feature_mode in ("curve", "curve_white"):
            cols = list(_curve_features(
                df["thetar"].to_numpy(), df["thetas"].to_numpy(),
                df["alpha_kpa"].to_numpy(), df["n"].to_numpy()).T)
        elif feature_mode == "vg":
            cols = [np.log10(df["alpha_kpa"]), np.log10(df["n"] - 1.0),
                    df["thetar"], df["thetas"]]
        else:
            raise ValueError(f"unknown feature_mode {feature_mode!r}")
        if use_depth:
            cols.append(np.log10(1.0 + df["depth_cm"].clip(lower=0)))
        feats = np.column_stack(cols)
        self.mean = feats.mean(axis=0)
        self.std = feats.std(axis=0)
        z = (feats - self.mean) / self.std
        # Curve features sampled at nearby heads are almost collinear (8 heads
        # carry only ~1.6 effective dimensions, mean |corr| 0.73), so plain
        # Euclidean distance over them collapses to "how wet is this soil" and
        # throws away curve shape. Whitening restores shape to equal footing.
        self._W = None
        if feature_mode == "curve_white":
            cov = np.cov(z, rowvar=False)
            ev, evec = np.linalg.eigh(cov)
            self._W = evec / np.sqrt(np.maximum(ev, 1e-8))
            z = z @ self._W
        self.z = z

    def set_excluded(self, layer_id):
        """Hide reference rows from all subsequent lookups. Accepts a single
        layer_id or a boolean mask over the reference rows (the latter is used
        to drop a whole profile at once). Used for large leave-one-out runs,
        where rebuilding the whole reference per target soil would be
        wasteful. Pass None to clear."""
        if layer_id is None:
            self._excluded = None
        elif isinstance(layer_id, np.ndarray) and layer_id.dtype == bool:
            self._excluded = layer_id
        else:
            self._excluded = (self.layer_id == layer_id)

    def set_excluded_profile(self, profile_id):
        """Hide every reference row belonging to one profile (all sibling
        horizons of a site), for grouped leave-one-out."""
        if self.profile_id is None:
            raise ValueError("reference has no profile_id column")
        self._excluded = (self.profile_id == profile_id)

    def neighbors(self, thetar, thetas, alpha, n, k, depth=None):
        if self.feature_mode in ("curve", "curve_white"):
            f = list(_curve_features(thetar, thetas, alpha, n).ravel())
        else:
            f = [np.log10(alpha), np.log10(n - 1.0), thetar, thetas]
        if self.use_depth:
            if depth is None:
                raise ValueError("this reference uses depth; pass depth=...")
            f.append(np.log10(1.0 + max(depth, 0.0)))
        f = (np.array(f) - self.mean) / self.std
        if self._W is not None:
            f = f @ self._W
        d = np.sqrt(((self.z - f) ** 2).sum(axis=1))
        excluded = getattr(self, "_excluded", None)
        if excluded is not None:
            d = np.where(excluded, np.inf, d)
        idx = np.argpartition(d, k)[:k]
        w = 1.0 / (d[idx] ** 2 + 1e-6)
        return idx, w / w.sum()


class TextureGBM:
    """Gradient-boosted classifier for the texture class only.

    Trained on the same four van Genuchten features the kNN uses. On matched
    folds it beats neighbour voting by +2.6 points in-distribution
    (p=0.019) and +5.4 points against an unseen laboratory (p=1.3e-07); see
    verify_gbm.py. It does not replace the kNN, which still supplies particle
    fractions, Ksat, the prediction intervals and the list of similar real
    soils -- a tree ensemble gives none of those.

    Training takes about 3 s on the 9,996-layer reference (early stopping
    settles around 50 iterations), so the model is fitted on demand rather
    than shipped as a pickle, which would tie the repository to one scikit-
    learn version.
    """

    def __init__(self, path=None, df=None, reference=DEFAULT_REFERENCE,
                 use_depth=False, random_state=0):
        try:
            from sklearn.ensemble import HistGradientBoostingClassifier
        except ImportError:
            raise SystemExit(
                "the hybrid model needs scikit-learn: pip install scikit-learn")
        import pandas as pd
        if df is None:
            df = (pd.read_csv(path) if path is not None
                  else load_reference_df(reference))
        if use_depth:
            df = df[df["depth_cm"].notna()]
        df = df[df["texture_class"].notna()]
        self.use_depth = use_depth
        self.clf = HistGradientBoostingClassifier(
            max_iter=400, learning_rate=0.06, max_leaf_nodes=31,
            l2_regularization=1.0, early_stopping=True,
            validation_fraction=0.15, random_state=random_state,
            # matches the uniform prior the kNN imposes at tau = 1
            class_weight="balanced")
        self.clf.fit(self._features(
            df["thetar"].to_numpy(), df["thetas"].to_numpy(),
            df["alpha_kpa"].to_numpy(), df["n"].to_numpy(),
            df["depth_cm"].to_numpy() if use_depth else None),
            df["texture_class"].to_numpy())
        self.classes_ = self.clf.classes_

    def _features(self, thetar, thetas, alpha, n, depth=None):
        cols = [np.log10(alpha), np.log10(np.asarray(n) - 1.0), thetar, thetas]
        if self.use_depth:
            cols.append(np.log10(1.0 + np.clip(np.asarray(depth, float), 0,
                                               None)))
        return np.column_stack(cols)

    def probabilities(self, thetar, thetas, alpha, n, depth=None):
        """Mean class probabilities over a set of Monte Carlo draws.

        Averaging predict_proba across the draws propagates the vG fit
        uncertainty exactly as the kNN's per-draw voting does.
        """
        d = (np.full(np.shape(thetar), depth if depth is not None else np.nan)
             if self.use_depth else None)
        p = self.clf.predict_proba(self._features(thetar, thetas, alpha, n, d))
        return dict(zip(self.classes_, p.mean(axis=0)))


def estimate(h, theta, ref=None, n_mc=300, k=30, seed=0, depth=None, clf=None):
    """Full pipeline: fit vG, then kNN inference with MC uncertainty.

    Returns a dict with fitted parameters, class probabilities, particle
    fractions and Ks with 5-95 % ranges.
    """
    h = np.asarray(h, dtype=float)
    theta = np.asarray(theta, dtype=float)
    if np.any(theta > 1.0):          # auto-detect percent input
        theta = theta / 100.0
    if ref is None:
        ref = GshpReference()

    popt, pcov = fit_vg(h, theta)
    thetar, thetas, la, ln1 = popt
    alpha, n = 10.0 ** la, 1.0 + 10.0 ** ln1
    perr = np.sqrt(np.diag(pcov))

    rng = np.random.default_rng(seed)
    draws = rng.multivariate_normal(popt, pcov, size=n_mc)
    draws = np.clip(draws, [0.0, 0.05, -4.0, np.log10(0.01)],
                    [0.5, 1.0, 1.5, np.log10(10.0)])

    votes = {c: 0.0 for c in USDA_CLASSES}
    frac_vals, frac_w = [], []
    ks_vals, ks_w = [], []
    ks_neighbor_counts = []
    for tr, ts, la_i, ln1_i in draws:
        idx, w = ref.neighbors(tr, ts, 10.0 ** la_i, 1.0 + 10.0 ** ln1_i, k,
                               depth=depth)
        w = w * ref.class_weight[idx] * ref.row_weight[idx]
        w = w / w.sum()
        ks_neighbor_counts.append(int(np.isfinite(ref.ksat[idx]).sum()))
        for i, wi in zip(idx, w):
            votes[ref.classes[i]] += wi
            frac_vals.append(ref.fractions[i]); frac_w.append(wi)
            if np.isfinite(ref.ksat[i]):
                ks_vals.append(ref.ksat[i]); ks_w.append(wi)

    total = sum(votes.values())
    knn_probs = {c: v / total for c, v in sorted(votes.items(),
                 key=lambda kv: -kv[1]) if v > 0}
    if clf is None:
        probs, source = knn_probs, "knn"
    else:
        # Hybrid: the class comes from the GBM, everything else below still
        # comes from the kNN neighbourhood computed above.
        p = clf.probabilities(draws[:, 0], draws[:, 1], 10.0 ** draws[:, 2],
                              1.0 + 10.0 ** draws[:, 3], depth=depth)
        probs = {c: v for c, v in sorted(p.items(), key=lambda kv: -kv[1])
                 if v > 0}
        source = "gbm"

    frac_vals = np.array(frac_vals); frac_w = np.array(frac_w)
    frac_mean = (frac_vals * frac_w[:, None]).sum(axis=0) / frac_w.sum()
    frac_lo = [_wpercentile(frac_vals[:, j], frac_w, 5) for j in range(3)]
    frac_hi = [_wpercentile(frac_vals[:, j], frac_w, 95) for j in range(3)]

    ks_vals = np.array(ks_vals); ks_w = np.array(ks_w)
    ks = {"n_neighbors_with_ksat": float(np.mean(ks_neighbor_counts)), "k": k}
    if len(ks_vals):
        log_ks = np.log10(ks_vals)
        ks.update(median_cmh=10.0 ** _wpercentile(log_ks, ks_w, 50),
                  p5_cmh=10.0 ** _wpercentile(log_ks, ks_w, 5),
                  p95_cmh=10.0 ** _wpercentile(log_ks, ks_w, 95))
    else:
        # No neighbor carried a measured Ksat -- report no estimate rather
        # than inventing one.
        ks.update(median_cmh=float("nan"), p5_cmh=float("nan"),
                  p95_cmh=float("nan"))

    return {
        "vg_fit": {
            "thetar": thetar, "thetas": thetas,
            "alpha_kpa": alpha, "n": n, "m": 1.0 - 1.0 / n,
            "se_thetar": perr[0], "se_thetas": perr[1],
            "se_log10_alpha": perr[2], "se_log10_n_minus_1": perr[3],
            "rmse": float(np.sqrt(np.mean(
                (vg_theta(h, thetar, thetas, alpha, n) - theta) ** 2))),
        },
        "texture_class": next(iter(probs)),
        "class_probabilities": probs,
        "class_source": source,
        "knn_class_probabilities": knn_probs,
        "fractions": {
            "sand": frac_mean[0], "silt": frac_mean[1], "clay": frac_mean[2],
            "p5": dict(zip(("sand", "silt", "clay"), frac_lo)),
            "p95": dict(zip(("sand", "silt", "clay"), frac_hi)),
            "class_of_mean": usda_class(*frac_mean),
        },
        "ksat": ks,
    }


def _wpercentile(x, w, q):
    order = np.argsort(x)
    cw = np.cumsum(w[order])
    return float(np.interp(q / 100.0 * cw[-1], cw, x[order]))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _read_data(path):
    rows = []
    with open(path, encoding="utf-8-sig") as f:
        for line in f:
            parts = line.replace(",", " ").split()
            if len(parts) >= 2:
                try:
                    rows.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    continue  # header line
    if not rows:
        sys.exit(f"no numeric (h, theta) rows found in {path}")
    return np.array(rows).T


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("datafile", help="two-column file: h_kPa, theta (fraction or %%)")
    ap.add_argument("--json", metavar="FILE", help="also write full results as JSON")
    ap.add_argument("--depth", type=float, metavar="CM",
                    help="sample mid-depth in cm; used as an extra matching "
                         "feature. Raises exact-class accuracy by ~3.5 points "
                         "on average, but only where the reference is dense "
                         "for the soil's region -- it LOWERS accuracy by ~7 "
                         "points for European soils, which are sparsely "
                         "represented. See README.")
    ap.add_argument("--model", choices=["knn", "hybrid"], default="knn",
                    help="knn (default) predicts the texture class by "
                         "neighbour voting. hybrid predicts it with a "
                         "gradient-boosted classifier, which is significantly "
                         "more accurate (+2.6 points in-distribution "
                         "p=0.019, +5.4 points against an unseen laboratory "
                         "p=1.3e-07) and is the better choice for a soil from "
                         "a new source. Fractions, Ksat, all intervals and "
                         "the neighbour list still come from the kNN either "
                         "way. Needs scikit-learn and ~3 s to train.")
    ap.add_argument("--tau", type=float, default=1.0, metavar="T",
                    help="class-prior exponent for neighbour votes, weight = "
                         "(1/class_count)^T. 1.0 (default) imposes a uniform "
                         "prior over the 12 classes and maximises recall for "
                         "rare classes, but over-predicts them: silt is "
                         "returned 2.6x more often than it occurs "
                         "(precision 20%%). 0.75 gives the best macro-F1 and "
                         "0.5 roughly calibrates silt, at almost no cost in "
                         "recall. 0 disables the correction. See README.")
    ap.add_argument("--reference", choices=sorted(REFERENCE_SETS),
                    default=DEFAULT_REFERENCE,
                    help="reference table. merged (default) = GSHP 9,996 "
                         "layers + NCSS/KSSL 2,530; gshp = GSHP only; "
                         "kssl = KSSL only (2,530 layers, and NO measured "
                         "Ksat at all, so Ks cannot be estimated); all adds "
                         "the 588 UNSODA 2.0 soils on top of merged. NOTE: "
                         "before 2026-08 the default was gshp and 'merged' "
                         "meant GSHP+UNSODA -- scripted callers should pass "
                         "--reference explicitly.")
    args = ap.parse_args()

    h, theta = _read_data(args.datafile)

    df = load_reference_df(args.reference)
    ref = GshpReference(df=df, use_depth=args.depth is not None, tau=args.tau)
    if df["ksat_cmh"].notna().sum() == 0:
        print(f"note: the '{args.reference}' reference carries no measured "
              f"Ksat, so no Ks estimate can be made.", file=sys.stderr)

    clf = None
    if args.model == "hybrid":
        # Trained on the same table the kNN uses, so the two halves of the
        # hybrid always see the same reference.
        clf = TextureGBM(df=df, use_depth=args.depth is not None)

    res = estimate(h, theta, ref=ref, depth=args.depth, clf=clf)

    vg = res["vg_fit"]
    print(f"van Genuchten fit (m = 1-1/n, alpha in kPa^-1):")
    print(f"  thetar = {vg['thetar']:.4f}  thetas = {vg['thetas']:.4f}  "
          f"alpha = {vg['alpha_kpa']:.4f}  n = {vg['n']:.3f}  m = {vg['m']:.3f}"
          f"  (RMSE {vg['rmse']:.4f})")
    src = ("gradient-boosted classifier" if res["class_source"] == "gbm"
           else "kNN neighbour vote")
    print(f"\nPredicted USDA texture class: {res['texture_class'].upper()}"
          f"   (from the {src})")
    for c, p in list(res["class_probabilities"].items())[:5]:
        print(f"  {c:<16s} {100 * p:5.1f} %")
    if res["class_source"] == "gbm":
        # The kNN's own answer is a free second opinion: agreement is a
        # meaningful confidence signal, disagreement a warning.
        knn_top = next(iter(res["knn_class_probabilities"]))
        agree = "agrees" if knn_top == res["texture_class"] else "DISAGREES"
        print(f"  kNN second opinion {agree}: {knn_top} "
              f"({100 * res['knn_class_probabilities'][knn_top]:.1f} %)")
    fr = res["fractions"]
    print(f"\nParticle fractions, mean [5-95 %]  (from the kNN; "
          f"mean plots as: {fr['class_of_mean']}):")
    for c in ("sand", "silt", "clay"):
        print(f"  {c:<5s} {fr[c]:5.1f} %  [{fr['p5'][c]:5.1f} - {fr['p95'][c]:5.1f}]")
    ks = res["ksat"]
    print(f"\nKs = {ks['median_cmh']:.3g} cm/h  "
          f"[5-95 %: {ks['p5_cmh']:.3g} - {ks['p95_cmh']:.3g}]  "
          f"(from ~{ks['n_neighbors_with_ksat']:.0f} of {ks['k']} "
          f"neighbors with measured Ksat)")

    if args.json:
        with open(args.json, "w") as f:
            json.dump(res, f, indent=2, default=float)
        print(f"\nfull results written to {args.json}")


if __name__ == "__main__":
    main()
