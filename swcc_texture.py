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
    # Laikipia, Arizona, the Canary Islands and the Yellow River Basin joined
    # the default in 2026-09, after EU-HYDI. EU-HYDI and Laikipia are restricted, so where their
    # tables have not been built the default quietly becomes "public" -- see
    # RESTRICTED_TABLES.
    "merged": ["gshp_reference.csv", "kssl_reference.csv",
               "hohenbrink_reference.csv", "babaeian_zanjanrood_reference.csv",
               "euhydi_reference.csv", "willard_reference.csv",
               "babaeian_az_reference.csv", "armas_reference.csv",
               "tong_reference.csv", "local_reference.csv"],
    # The distributed tables: what a fresh clone has.
    "public": ["gshp_reference.csv", "kssl_reference.csv",
               "hohenbrink_reference.csv", "babaeian_zanjanrood_reference.csv",
               "babaeian_az_reference.csv", "armas_reference.csv",
               "tong_reference.csv"],
    "gshp": ["gshp_reference.csv"],
    "kssl": ["kssl_reference.csv"],
    "hohenbrink": ["hohenbrink_reference.csv"],
    "euhydi": ["euhydi_reference.csv"],
    # Laikipia, Kenya (Willard et al., unpublished; prepare_willard.py); in
    # the default since 2026-09 after verify_external.py, without its field Ks.
    "willard": ["willard_reference.csv"],
    # Zanjanrood, Iran (Babaeian et al. 2015; prepare_babaeian_zanjanrood.py);
    # in the default since 2026-09 after verify_external.py.
    "babaeian_zanjanrood": ["babaeian_zanjanrood_reference.csv"],
    # Arizona, USA (Babaeian, unpublished; prepare_babaeian_az.py); in the
    # default since 2026-09 after verify_external.py.
    "babaeian_az": ["babaeian_az_reference.csv"],
    # Canary Islands andic soils (Armas Espinel 2013; prepare_armas.py); in
    # the default since 2026-09 after verify_external.py and verify_andic.py.
    "armas": ["armas_reference.csv"],
    # Yellow River Basin, China (Tong et al. 2024; prepare_tong.py); in the
    # default since 2026-09: read poorly as a new source (centrifuge curves,
    # laser-diffraction texture) but no harm to the rest and +6 pp for silt
    # loam when added (verify_external.py).
    "tong": ["tong_reference.csv"],
    "all": ["gshp_reference.csv", "kssl_reference.csv",
            "hohenbrink_reference.csv", "babaeian_zanjanrood_reference.csv",
            "euhydi_reference.csv", "willard_reference.csv",
            "babaeian_az_reference.csv", "armas_reference.csv",
            "tong_reference.csv", "unsoda_reference.csv", "sdb_reference.csv",
            "local_reference.csv"],
}
DEFAULT_REFERENCE = "merged"


# Tables built locally from restricted data and never distributed. A named set
# that includes one still loads without it -- from its other tables, with a
# note on stderr -- so a fresh clone of the public repository runs. A set made
# only of restricted tables cannot load at all.
RESTRICTED_TABLES = {"euhydi_reference.csv", "willard_reference.csv"}
_noted_missing = set()

# A user's own lab-verified curves, written by add_local_data.py. Joins the
# default reference whenever it exists; its absence is normal and silent.
OPTIONAL_TABLES = {"local_reference.csv"}

# Volcanic parent material and andic properties per layer (prepare_volcanic.py),
# joined on layer_id as the columns volcanic and andic when the tables exist.
VOLCANIC_TABLES = ("volcanic_flags.csv", "volcanic_flags_restricted.csv")

# Andic soil properties as a user input (--andic yes/no). The classifier takes
# it as a feature (reference: yes/likely 1, no 0, unknown missing), and the
# neighbour search adds ANDIC_LAMBDA standard units of distance between soils
# that differ in it (reference unknown counts as not andic) -- for the class
# vote and the fractions only. Ks neighbours ignore it: when tested, the only
# andic soils with a measured Ks came from one source (the Canary Islands
# have since added a second). See verify_volcanic.py and verify_andic_knn.py.
ANDIC_CODE = {"yes": 1.0, "likely": 1.0, "no": 0.0}
ANDIC_LAMBDA = 1.0

# How the wet end of a sample's retention curve was measured: on an intact
# core or on repacked material (see the prepare_*.py scripts). Encoded for the
# GBM; "unknown" and missing become NaN, which the GBM handles natively.
SAMPLE_TYPE_CODE = {"undisturbed": 0.0, "disturbed": 1.0}

# Local covariates a user may supply: sample mid-depth (cm) and bulk density
# (g/cm3). Reference column names; see TextureGBM(covariates=...).
COVARIATES = ("depth_cm", "bd")


def load_reference_df(name=DEFAULT_REFERENCE):
    """Concatenate the CSVs making up a named reference set."""
    import pandas as pd
    if name not in REFERENCE_SETS:
        raise ValueError(f"unknown reference set {name!r}; "
                         f"choose from {sorted(REFERENCE_SETS)}")
    files = REFERENCE_SETS[name]
    have = [f for f in files if os.path.exists(os.path.join(DATA_DIR, f))]
    missing = [f for f in files if f not in have and f not in OPTIONAL_TABLES]
    if not have or any(f not in RESTRICTED_TABLES for f in missing):
        raise FileNotFoundError(
            f"reference set {name!r} needs {', '.join(missing)}, which is not "
            f"in this checkout. Build it with the matching prepare_*.py; "
            f"restricted tables are not distributed.")
    for f in missing:
        if f not in _noted_missing:
            _noted_missing.add(f)
            print(f"note: {f} is built locally from restricted data and is "
                  f"not in this checkout, so the {name!r} reference uses the "
                  f"public tables only.", file=sys.stderr)
    df = pd.concat([pd.read_csv(os.path.join(DATA_DIR, f)) for f in have],
                   ignore_index=True)
    flags = [pd.read_csv(os.path.join(DATA_DIR, f),
                         usecols=["layer_id", "volcanic", "andic"])
             for f in VOLCANIC_TABLES if os.path.exists(os.path.join(DATA_DIR, f))]
    if flags:
        f = pd.concat(flags, ignore_index=True)
        f["_k"] = f.layer_id.astype(str)
        df["_k"] = df.layer_id.astype(str)
        df = df.merge(f.drop(columns="layer_id").drop_duplicates("_k"),
                      on="_k", how="left").drop(columns="_k")
        df["andic_code"] = df["andic"].map(ANDIC_CODE)
    return df

USDA_CLASSES = [
    "sand", "loamy sand", "sandy loam", "loam", "silt", "silt loam",
    "sandy clay loam", "clay loam", "silty clay loam", "sandy clay",
    "silty clay", "clay",
]

# ---------------------------------------------------------------------------
# van Genuchten model and fitting
# ---------------------------------------------------------------------------

def vg_theta(h, thetar, thetas, alpha, n, m=None):
    """van Genuchten retention, h in kPa >= 0.

    m defaults to the Mualem constraint 1 - 1/n. Passing m explicitly gives
    the unconstrained (five-parameter) form, in which n and m independently
    control the two inflection regions of the curve -- the dry end, where
    texture information lives, is governed mostly by m.
    """
    if m is None:
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


def _vg_transformed_m(h, thetar, thetas, la, ln1, m):
    return vg_theta(h, thetar, thetas, 10.0 ** la, 1.0 + 10.0 ** ln1, m)


# Bounds on the free m. Away from 0 and 1 because both ends are degenerate:
# m -> 0 flattens the curve entirely, m -> 1 with large n makes the dry limb
# vertical, and curve_fit wanders into both if allowed.
M_BOUNDS = (0.02, 0.98)


def fit_vg(h, theta, free_m=False):
    """Fit vG parameters to (h [kPa], theta [m3/m3]) data.

    Fits in transformed space p = [thetar, thetas, log10(alpha), log10(n-1)]
    for stability, with m tied to n by the Mualem constraint. With
    free_m=True a fifth parameter m is fitted independently and appended to
    p, which needs at least 6 points. Returns (popt, pcov) in that space.
    """
    h = np.asarray(h, dtype=float)
    theta = np.asarray(theta, dtype=float)
    npar = 5 if free_m else 4
    if len(h) < npar + 1:
        raise ValueError(f"need at least {npar + 1} (h, theta) points to fit "
                         f"{npar} parameters")

    tmin, tmax = theta.min(), theta.max()
    hpos = h[h > 0]
    p0 = [max(0.5 * tmin, 1e-3), tmax, np.log10(1.0 / np.median(hpos)), np.log10(0.5)]
    lb = [0.0, 0.5 * tmax, -4.0, np.log10(0.01)]
    ub = [tmin + 1e-9, 1.0, 1.5, np.log10(10.0)]
    f = _vg_transformed
    if free_m:
        # Start from the Mualem value implied by p0's n, so the free fit
        # begins at the constrained solution and only departs if the data
        # ask it to.
        p0 = p0 + [1.0 - 1.0 / (1.0 + 10.0 ** p0[3])]
        lb = lb + [M_BOUNDS[0]]
        ub = ub + [M_BOUNDS[1]]
        f = _vg_transformed_m
    p0 = np.clip(p0, lb, ub)

    popt, pcov = curve_fit(f, h, theta, p0=p0, bounds=(lb, ub), maxfev=20000)
    return popt, pcov


# ---------------------------------------------------------------------------
# USDA texture triangle
# ---------------------------------------------------------------------------

def usda_class(sand, silt, clay):
    """USDA texture class from fractions in % (must sum to ~100).

    Boundaries follow the NRCS Soil Texture Calculator: a class's lower clay
    bound is inclusive, so exactly 20 % clay is sandy clay loam (not sandy
    loam) and exactly 27 % is clay loam (not loam).
    """
    if silt + 1.5 * clay < 15:
        return "sand"
    if silt + 2.0 * clay < 30:
        return "loamy sand"
    if (7 <= clay < 20 and sand > 52 and silt + 2 * clay >= 30) or \
       (clay < 7 and silt < 50 and silt + 2 * clay >= 30):
        return "sandy loam"
    if 7 <= clay < 27 and 28 <= silt < 50 and sand <= 52:
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
                 use_depth=False, use_om=False, use_bd=False, use_m=False,
                 tau=1.0, feature_mode="vg",
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
        self.use_om = use_om
        self.use_bd = use_bd
        # use_m switches to the unconstrained five-parameter vG coordinates
        # (see free_m.py): freeing m from n gives the dry limb -- where
        # texture information sits -- its own axis.
        self.use_m = use_m
        if use_depth:
            df = df[df["depth_cm"].notna()].reset_index(drop=True)
        if use_om:
            # Organic carbon is present for only ~32 % of the merged
            # reference, so switching it on discards most of the rows. Whether
            # the covariate repays that is measured in verify_om.py.
            df = df[df["oc"].notna()].reset_index(drop=True)
        if use_bd:
            # Bulk density (oven-dry, g/cm3) as an extra matching dimension,
            # used when the user supplies it. Rows without a plausible value
            # (~0.5 %) are dropped. It lowers the Ks error slightly when the
            # user's own source is in the reference (0.42 -> 0.39 dex,
            # p=0.04); see verify_covariates.py.
            df = df[df["bd"].between(0.1, 2.3)].reset_index(drop=True)
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
        self.sample_type = (df["sample_type"].fillna("unknown").to_numpy()
                            if "sample_type" in df.columns else None)
        self.andic = (df["andic"].isin(["yes", "likely"]).to_numpy(float)
                      if "andic" in df.columns else None)
        if feature_mode in ("curve", "curve_white"):
            cols = list(_curve_features(
                df["thetar"].to_numpy(), df["thetas"].to_numpy(),
                df["alpha_kpa"].to_numpy(), df["n"].to_numpy()).T)
        elif feature_mode == "vg":
            if use_m:
                cols = [np.log10(df["alpha_kpa_m"]), np.log10(df["n_m"] - 1.0),
                        df["thetar_m"], df["thetas_m"], df["m"]]
            else:
                cols = [np.log10(df["alpha_kpa"]), np.log10(df["n"] - 1.0),
                        df["thetar"], df["thetas"]]
        else:
            raise ValueError(f"unknown feature_mode {feature_mode!r}")
        if use_depth:
            cols.append(np.log10(1.0 + df["depth_cm"].clip(lower=0)))
        if use_om:
            cols.append(np.log10(1.0 + df["oc"].clip(lower=0)))
        if use_bd:
            cols.append(np.log10(1.0 + df["bd"]))
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

    def neighbors(self, thetar, thetas, alpha, n, k, depth=None, om=None,
                  m=None, sample_type=None, bd=None, andic=None):
        """k nearest reference rows and their normalised 1/d^2 weights.
        With sample_type, only rows of that type are eligible. With andic
        (1.0 or 0.0), rows that differ from it are ANDIC_LAMBDA standard
        units further away."""
        if self.feature_mode in ("curve", "curve_white"):
            f = list(_curve_features(thetar, thetas, alpha, n).ravel())
        elif self.use_m:
            if m is None:
                raise ValueError("this reference uses free m; pass m=...")
            f = [np.log10(alpha), np.log10(n - 1.0), thetar, thetas, m]
        else:
            f = [np.log10(alpha), np.log10(n - 1.0), thetar, thetas]
        if self.use_depth:
            if depth is None:
                raise ValueError("this reference uses depth; pass depth=...")
            f.append(np.log10(1.0 + max(depth, 0.0)))
        if self.use_om:
            if om is None:
                raise ValueError("this reference uses organic carbon; pass om=")
            f.append(np.log10(1.0 + max(om, 0.0)))
        if self.use_bd:
            if bd is None:
                raise ValueError("this reference uses bulk density; pass bd=")
            f.append(np.log10(1.0 + bd))
        f = (np.array(f) - self.mean) / self.std
        if self._W is not None:
            f = f @ self._W
        d2 = ((self.z - f) ** 2).sum(axis=1)
        if andic is not None:
            if self.andic is None:
                raise ValueError("this reference has no andic flags; run "
                                 "prepare_volcanic.py")
            d2 = d2 + (ANDIC_LAMBDA * (self.andic - andic)) ** 2
        d = np.sqrt(d2)
        excluded = getattr(self, "_excluded", None)
        if excluded is not None:
            d = np.where(excluded, np.inf, d)
        if sample_type is not None:
            d = np.where(self.sample_type == sample_type, d, np.inf)
        idx = np.argpartition(d, k)[:k]
        w = 1.0 / (d[idx] ** 2 + 1e-6)
        return idx, w / w.sum()


class TextureGBM:
    """Gradient-boosted classifier for the texture class only.

    Trained on the same four van Genuchten features the kNN uses, plus the
    sample type. On the 20,052-layer reference it is level with neighbour
    voting on exact class, in-distribution (-0.5 points, p=0.72) and against
    an unseen laboratory (-0.2, p=0.88), and better on the texture group for
    an unseen laboratory (+2.3 points); see verify_hybrid.py. It
    does not replace the kNN, which still supplies particle
    fractions, Ksat, the prediction intervals and the list of similar real
    soils -- a tree ensemble gives none of those.

    Training takes about 3 s on the 9,996-layer reference (early stopping
    settles around 50 iterations), so the model is fitted on demand rather
    than shipped as a pickle, which would tie the repository to one scikit-
    learn version.
    """

    def __init__(self, path=None, df=None, reference=DEFAULT_REFERENCE,
                 use_depth=False, use_m=False, use_sample_type=True,
                 covariates=(), random_state=0):
        try:
            from sklearn.ensemble import HistGradientBoostingClassifier
        except ImportError:
            raise ImportError(
                "the hybrid model needs scikit-learn: pip install scikit-learn")
        import pandas as pd
        if df is None:
            df = (pd.read_csv(path) if path is not None
                  else load_reference_df(reference))
        if use_depth:
            df = df[df["depth_cm"].notna()]
        df = df[df["texture_class"].notna()]
        self.use_depth = use_depth
        self.use_m = use_m
        # sample_type as a fifth feature. For undisturbed soils it changes
        # little (+0.0 points, p=1.00, own laboratory in the reference; -0.9,
        # p=0.58, unseen laboratory); see verify_sample_type.py.
        self.use_sample_type = use_sample_type
        # Optional user covariates, a subset of COVARIATES, appended as extra
        # features. Reference rows lacking one keep it as NaN. Depth and bulk
        # density together add ~4 points when the user's own data source is
        # in the reference and change little (-1 to +2.5 points, depending on
        # the targets) for a source it has never seen; see
        # verify_covariates.py.
        self.covariates = tuple(covariates)
        self.clf = HistGradientBoostingClassifier(
            max_iter=400, learning_rate=0.06, max_leaf_nodes=31,
            l2_regularization=1.0, early_stopping=True,
            validation_fraction=0.15, random_state=random_state,
            # matches the uniform prior the kNN imposes at tau = 1
            class_weight="balanced")
        cols = ("thetar_m", "thetas_m", "alpha_kpa_m", "n_m") if use_m \
            else ("thetar", "thetas", "alpha_kpa", "n")
        stype = (df["sample_type"].map(SAMPLE_TYPE_CODE).to_numpy(float)
                 if use_sample_type and "sample_type" in df.columns
                 else np.full(len(df), np.nan))
        covs = [df[c].to_numpy(float) if c in df.columns
                else np.full(len(df), np.nan) for c in self.covariates]
        self.clf.fit(self._features(
            *(df[c].to_numpy() for c in cols),
            df["depth_cm"].to_numpy() if use_depth else None,
            df["m"].to_numpy() if use_m else None, stype, covs),
            df["texture_class"].to_numpy())
        self.classes_ = self.clf.classes_

    def _features(self, thetar, thetas, alpha, n, depth=None, m=None,
                  stype=None, covs=()):
        cols = [np.log10(alpha), np.log10(np.asarray(n) - 1.0), thetar, thetas]
        if self.use_depth:
            cols.append(np.log10(1.0 + np.clip(np.asarray(depth, float), 0,
                                               None)))
        if self.use_m:
            cols.append(np.asarray(m, dtype=float))
        if self.use_sample_type:
            cols.append(np.asarray(stype, dtype=float))
        cols.extend(np.asarray(c, dtype=float) for c in covs)
        return np.column_stack(cols)

    def probabilities(self, thetar, thetas, alpha, n, depth=None, m=None,
                      sample_type=None, bulk_density=None, andic=None):
        """Mean class probabilities over a set of Monte Carlo draws.

        Averaging predict_proba across the draws propagates the vG fit
        uncertainty exactly as the kNN's per-draw voting does.
        """
        d = (np.full(np.shape(thetar), depth if depth is not None else np.nan)
             if self.use_depth else None)
        st_code = np.full(np.shape(thetar),
                          SAMPLE_TYPE_CODE.get(sample_type, np.nan))
        given = {"depth_cm": depth, "bd": bulk_density, "andic_code": andic}
        covs = [np.full(np.shape(thetar), np.nan if given[c] is None
                        else given[c], dtype=float) for c in self.covariates]
        p = self.clf.predict_proba(
            self._features(thetar, thetas, alpha, n, d, m, st_code, covs))
        return dict(zip(self.classes_, p.mean(axis=0)))


def estimate(h, theta, ref=None, n_mc=300, k=30, seed=0, depth=None,
             om=None, clf=None, sample_type=None, bulk_density=None,
             andic=None):
    """Full pipeline: fit vG, then kNN inference with MC uncertainty.

    Returns a dict with fitted parameters, class probabilities, particle
    fractions and Ks with 5-95 % ranges.

    sample_type ("undisturbed" or "disturbed") is how the target's curve was
    measured. The GBM receives it as a feature, and for "undisturbed" Ks comes
    only from undisturbed reference soils; the kNN class vote and the particle
    fractions still use the whole reference. None or "unknown" uses neither.

    depth (cm) and bulk_density (g/cm3) reach the classifier when it was
    trained with them (TextureGBM(covariates=...)); bulk_density also joins
    the neighbour search when the reference was built with use_bd=True.

    andic ("yes"/"no", or None for not known) reaches the classifier when it
    was trained with covariates=["andic_code"], and the neighbour search for
    the class vote and the fractions; Ks neighbours ignore it.
    """
    andic_code = (None if andic is None else
                  ANDIC_CODE[andic] if isinstance(andic, str) else float(bool(andic)))
    h = np.asarray(h, dtype=float)
    theta = np.asarray(theta, dtype=float)
    if np.any(theta > 1.0):          # auto-detect percent input
        theta = theta / 100.0
    if ref is None:
        ref = GshpReference()

    # The reference decides the coordinate system: a use_m reference is
    # indexed on the unconstrained five-parameter fit, so the target has to
    # be fitted the same way.
    free_m = bool(getattr(ref, "use_m", False))
    popt, pcov = fit_vg(h, theta, free_m=free_m)
    thetar, thetas, la, ln1 = popt[:4]
    alpha, n = 10.0 ** la, 1.0 + 10.0 ** ln1
    m_fit = popt[4] if free_m else 1.0 - 1.0 / n
    perr = np.sqrt(np.diag(pcov))

    rng = np.random.default_rng(seed)
    draws = rng.multivariate_normal(popt, pcov, size=n_mc)
    lo = [0.0, 0.05, -4.0, np.log10(0.01)]
    hi = [0.5, 1.0, 1.5, np.log10(10.0)]
    if free_m:
        lo, hi = lo + [M_BOUNDS[0]], hi + [M_BOUNDS[1]]
    draws = np.clip(draws, lo, hi)

    # An undisturbed sample takes Ks only from undisturbed neighbours: near
    # saturation an intact core keeps its macropores and a repacked one does
    # not. Disturbed samples use every soil, because the reference holds
    # almost no disturbed layers with a measured Ks (45 of 11,017). Matching
    # the class and fractions on type as well was tested and cost accuracy
    # against an unseen laboratory, so only Ks is restricted.
    ks_type = ("undisturbed" if sample_type == "undisturbed"
               and getattr(ref, "sample_type", None) is not None else None)

    votes = {c: 0.0 for c in USDA_CLASSES}
    frac_vals, frac_w = [], []
    ks_vals, ks_w = [], []
    ks_neighbor_counts = []
    for row in draws:
        tr, ts, la_i, ln1_i = row[:4]
        idx, w = ref.neighbors(tr, ts, 10.0 ** la_i, 1.0 + 10.0 ** ln1_i, k,
                               depth=depth, om=om, bd=bulk_density,
                               m=row[4] if free_m else None, andic=andic_code)
        w = w * ref.class_weight[idx] * ref.row_weight[idx]
        w = w / w.sum()
        for i, wi in zip(idx, w):
            votes[ref.classes[i]] += wi
            frac_vals.append(ref.fractions[i]); frac_w.append(wi)
        if ks_type is not None or andic_code is not None:
            idx, w = ref.neighbors(tr, ts, 10.0 ** la_i, 1.0 + 10.0 ** ln1_i,
                                   k, depth=depth, om=om, bd=bulk_density,
                                   m=row[4] if free_m else None,
                                   sample_type=ks_type)
            w = w * ref.class_weight[idx] * ref.row_weight[idx]
            w = w / w.sum()
        ks_neighbor_counts.append(int(np.isfinite(ref.ksat[idx]).sum()))
        for i, wi in zip(idx, w):
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
                              1.0 + 10.0 ** draws[:, 3], depth=depth,
                              m=draws[:, 4] if free_m else None,
                              sample_type=sample_type,
                              bulk_density=bulk_density, andic=andic_code)
        probs = {c: v for c, v in sorted(p.items(), key=lambda kv: -kv[1])
                 if v > 0}
        source = "gbm"

    frac_vals = np.array(frac_vals); frac_w = np.array(frac_w)
    frac_mean = (frac_vals * frac_w[:, None]).sum(axis=0) / frac_w.sum()
    frac_lo = [_wpercentile(frac_vals[:, j], frac_w, 5) for j in range(3)]
    frac_hi = [_wpercentile(frac_vals[:, j], frac_w, 95) for j in range(3)]

    ks_vals = np.array(ks_vals); ks_w = np.array(ks_w)
    ks = {"n_neighbors_with_ksat": float(np.mean(ks_neighbor_counts)), "k": k,
          "sample_type": ks_type}
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
            "alpha_kpa": alpha, "n": n, "m": m_fit,
            "m_is_free": free_m,
            "se_thetar": perr[0], "se_thetas": perr[1],
            "se_log10_alpha": perr[2], "se_log10_n_minus_1": perr[3],
            "rmse": float(np.sqrt(np.mean(
                (vg_theta(h, thetar, thetas, alpha, n,
                          m_fit if free_m else None) - theta) ** 2))),
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
                    help="sample mid-depth in cm. Together with "
                         "--bulk-density it adds ~4 points of exact class "
                         "when your own verified data are in the reference "
                         "(add_local_data.py); for a data source the "
                         "reference has never seen the effect is small and "
                         "not stable (-1 to +2.5 points). Part of the gain "
                         "favours the classes your own data hold most. See "
                         "README.")
    ap.add_argument("--bulk-density", type=float, metavar="G_CM3",
                    help="bulk density in g/cm3 (oven-dry). The classifier "
                         "uses it as for --depth, and the neighbour search "
                         "also matches on it, which lowers the Ks error "
                         "slightly. See README.")
    ap.add_argument("--andic", choices=["yes", "no"],
                    help="whether the soil has andic properties (an Andisol "
                         "or Andosol, or andic by oxalate Al + 1/2 Fe and "
                         "phosphate retention). The classifier uses it and "
                         "the neighbour search matches on it for the class "
                         "vote and the fractions, not for Ks. For andic soils "
                         "from a source the reference has never seen it adds "
                         "~6-13 points of exact class and ~13-20 of texture "
                         "group, and brings the fractions 1.5 points closer; "
                         "other soils do not change. Leave it out when not "
                         "known. See README.")
    ap.add_argument("--model", choices=["knn", "hybrid"], default=None,
                    help="hybrid (default) predicts the texture class with a "
                         "gradient-boosted classifier, which is level with "
                         "neighbour voting in-distribution (+0.6 points) and "
                         "better against an unseen laboratory (+1.2 points "
                         "exact class, +3.2 points texture group). "
                         "knn predicts it by "
                         "neighbour voting instead, which needs no "
                         "scikit-learn, runs in ~0.6 s rather than ~5 s, and "
                         "is bit-reproducible across environments. Fractions, "
                         "Ksat, all intervals and the neighbour list come from "
                         "the kNN either way; the hybrid also prints the kNN's "
                         "class as a second opinion. If scikit-learn is "
                         "missing the default falls back to knn.")
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
                         "layers + NCSS/KSSL 2,530 + Hohenbrink 560 + "
                         "Zanjanrood (Babaeian) 169 + EU-HYDI 6,797 + "
                         "Laikipia 86 + Arizona 21 + Canary Islands 66 + "
                         "Yellow River 1,030, 21,255 "
                         "in total. EU-HYDI (consortium-restricted) and Laikipia "
                         "(unpublished) are built locally with "
                         "prepare_euhydi.py and prepare_willard.py and never "
                         "distributed, so where they are absent merged uses "
                         "the other seven (14,372) and says so. public = those "
                         "seven distributed tables; gshp, kssl, "
                         "hohenbrink, babaeian_zanjanrood, euhydi, willard, "
                         "babaeian_az, armas and tong select "
                         "one source "
                         "(kssl has NO measured Ksat, so Ks cannot be "
                         "estimated from it); all adds UNSODA 2.0 and sDB. "
                         "NOTE: what merged contains has changed over time -- "
                         "scripted callers should pass --reference explicitly.")
    ap.add_argument("--sample-type", default="undisturbed",
                    choices=["undisturbed", "disturbed", "unknown"],
                    help="how the curve was measured. undisturbed (default) "
                         "suits in-situ sensors and intact cores -- Ks should "
                         "be measured on intact soil, where the large pores "
                         "survive; disturbed is for repacked or sieved "
                         "samples; unknown gives no type. The classifier is "
                         "told the type, and an undisturbed sample takes Ks "
                         "only from undisturbed reference soils (10,951 of the "
                         "11,017 with a measured Ks). The reference has almost "
                         "no disturbed Ks, so disturbed and unknown samples "
                         "take Ks from all soils. See README.")
    args = ap.parse_args()

    h, theta = _read_data(args.datafile)

    df = load_reference_df(args.reference)
    if args.andic is not None and "andic_code" not in df.columns:
        sys.exit("--andic needs the volcanic flags (data/volcanic_flags.csv); "
                 "build them with prepare_volcanic.py")
    ref = GshpReference(df=df, tau=args.tau,
                        use_bd=args.bulk_density is not None)
    if df["ksat_cmh"].notna().sum() == 0:
        print(f"note: the '{args.reference}' reference carries no measured "
              f"Ksat, so no Ks estimate can be made.", file=sys.stderr)

    clf = None
    model = args.model or "hybrid"
    if model == "hybrid":
        try:
            # Trained on the same table the kNN uses, so the two halves of the
            # hybrid always see the same reference.
            given = {"depth_cm": args.depth, "bd": args.bulk_density,
                     "andic_code": args.andic}
            clf = TextureGBM(df=df, covariates=[c for c in given
                                                if given[c] is not None])
        except ImportError as exc:
            if args.model is not None:
                raise SystemExit(str(exc))
            # The hybrid is only the default, not a request: degrade to the
            # kNN rather than refusing to run.
            print(f"note: {exc}; falling back to --model knn.", file=sys.stderr)

    if args.sample_type == "disturbed":
        n_dis = int(((ref.sample_type == "disturbed")
                     & np.isfinite(ref.ksat)).sum())
        print(f"note: the reference holds only {n_dis} disturbed layers with a "
              f"measured Ks, so Ks comes from all soils.", file=sys.stderr)
    if clf is None and args.depth is not None:
        print("note: --depth is used by the classifier only, so it is ignored "
              "with --model knn.", file=sys.stderr)
    res = estimate(h, theta, ref=ref, depth=args.depth, clf=clf,
                   sample_type=args.sample_type,
                   bulk_density=args.bulk_density, andic=args.andic)

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
    print(f"\nParticle fractions, mean [5-95 %]  (from the kNN"
          f"{'; matched on andic properties' if args.andic else ''}; "
          f"mean plots as: {fr['class_of_mean']}):")
    for c in ("sand", "silt", "clay"):
        print(f"  {c:<5s} {fr[c]:5.1f} %  [{fr['p5'][c]:5.1f} - {fr['p95'][c]:5.1f}]")
    ks = res["ksat"]
    print(f"\nKs = {ks['median_cmh']:.3g} cm/h  "
          f"[5-95 %: {ks['p5_cmh']:.3g} - {ks['p95_cmh']:.3g}]  "
          f"(from ~{ks['n_neighbors_with_ksat']:.0f} of {ks['k']} "
          f"{ks['sample_type'] + ' ' if ks['sample_type'] else ''}"
          f"neighbors with measured Ksat"
          f"{', matched on bulk density' if ref.use_bd else ''})")

    if args.json:
        with open(args.json, "w") as f:
            json.dump(res, f, indent=2, default=float)
        print(f"\nfull results written to {args.json}")


if __name__ == "__main__":
    main()
