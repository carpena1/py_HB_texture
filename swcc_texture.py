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

# Resolve the reference table relative to this script, so the tool works
# from any working directory.
REFERENCE_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "data", "gshp_reference.csv")

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
    def __init__(self, path=REFERENCE_CSV, df=None, use_depth=False):
        """Build the reference from the CSV at `path`, or from a preloaded
        DataFrame `df` (used for leave-one-out verification, where the target
        soil is dropped before constructing the reference).

        `use_depth=True` adds the sample mid-depth as a fifth, equally
        weighted feature (log10(1+depth_cm), standardized); reference rows
        without a depth are dropped in that case.
        """
        if df is None:
            import pandas as pd
            df = pd.read_csv(path)
        self.use_depth = use_depth
        if use_depth:
            df = df[df["depth_cm"].notna()].reset_index(drop=True)
        self.layer_id = df["layer_id"].to_numpy()
        self.profile_id = (df["profile_id"].to_numpy()
                           if "profile_id" in df.columns else None)
        self.classes = df["texture_class"].to_numpy()
        # Inverse class-frequency weights: votes assume a uniform prior over
        # the 12 USDA classes rather than GSHP's sampling distribution
        # (sand alone is ~39 % of the database).
        counts = df["texture_class"].value_counts()
        self.class_weight = (1.0 / counts)[df["texture_class"]].to_numpy()
        self.fractions = df[["sand", "silt", "clay"]].to_numpy()
        self.ksat = df["ksat_cmh"].to_numpy()
        cols = [np.log10(df["alpha_kpa"]), np.log10(df["n"] - 1.0),
                df["thetar"], df["thetas"]]
        if use_depth:
            cols.append(np.log10(1.0 + df["depth_cm"].clip(lower=0)))
        feats = np.column_stack(cols)
        self.mean = feats.mean(axis=0)
        self.std = feats.std(axis=0)
        self.z = (feats - self.mean) / self.std

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
        f = [np.log10(alpha), np.log10(n - 1.0), thetar, thetas]
        if self.use_depth:
            if depth is None:
                raise ValueError("this reference uses depth; pass depth=...")
            f.append(np.log10(1.0 + max(depth, 0.0)))
        f = np.array(f)
        d = np.sqrt(((self.z - (f - self.mean) / self.std) ** 2).sum(axis=1))
        excluded = getattr(self, "_excluded", None)
        if excluded is not None:
            d = np.where(excluded, np.inf, d)
        idx = np.argpartition(d, k)[:k]
        w = 1.0 / (d[idx] ** 2 + 1e-6)
        return idx, w / w.sum()


def estimate(h, theta, ref=None, n_mc=300, k=30, seed=0, depth=None):
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
        w = w * ref.class_weight[idx]
        w = w / w.sum()
        ks_neighbor_counts.append(int(np.isfinite(ref.ksat[idx]).sum()))
        for i, wi in zip(idx, w):
            votes[ref.classes[i]] += wi
            frac_vals.append(ref.fractions[i]); frac_w.append(wi)
            if np.isfinite(ref.ksat[i]):
                ks_vals.append(ref.ksat[i]); ks_w.append(wi)

    total = sum(votes.values())
    probs = {c: v / total for c, v in sorted(votes.items(),
             key=lambda kv: -kv[1]) if v > 0}

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
    ap.add_argument("--reference", choices=["gshp", "merged", "kssl", "all"],
                    default="gshp",
                    help="reference database: gshp (default); merged adds the "
                         "588 UNSODA 2.0 soils; kssl adds 2,508 NCSS/KSSL "
                         "layers; all adds both. Neither merge reaches "
                         "significance on its own (UNSODA +0.7 points "
                         "p=0.052, KSSL +1.2 points p=0.090), and KSSL "
                         "carries no Ksat. See README.")
    args = ap.parse_args()

    h, theta = _read_data(args.datafile)

    ref = None
    if args.reference != "gshp" or args.depth is not None:
        import pandas as pd
        here = os.path.dirname(os.path.abspath(__file__))
        df = pd.read_csv(REFERENCE_CSV)
        extra = {"merged": ["unsoda_reference.csv"],
                 "kssl": ["kssl_reference.csv"],
                 "all": ["unsoda_reference.csv", "kssl_reference.csv"]}
        for name in extra.get(args.reference, []):
            df = pd.concat([df, pd.read_csv(os.path.join(here, "data", name))],
                           ignore_index=True)
        ref = GshpReference(df=df, use_depth=args.depth is not None)

    res = estimate(h, theta, ref=ref, depth=args.depth)

    vg = res["vg_fit"]
    print(f"van Genuchten fit (m = 1-1/n, alpha in kPa^-1):")
    print(f"  thetar = {vg['thetar']:.4f}  thetas = {vg['thetas']:.4f}  "
          f"alpha = {vg['alpha_kpa']:.4f}  n = {vg['n']:.3f}  m = {vg['m']:.3f}"
          f"  (RMSE {vg['rmse']:.4f})")
    print(f"\nPredicted USDA texture class: {res['texture_class'].upper()}")
    for c, p in list(res["class_probabilities"].items())[:5]:
        print(f"  {c:<16s} {100 * p:5.1f} %")
    fr = res["fractions"]
    print(f"\nParticle fractions, mean [5-95 %]  "
          f"(mean plots as: {fr['class_of_mean']}):")
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
