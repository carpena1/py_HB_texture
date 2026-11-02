"""Shared constants and statistics for the verification scripts.

Collected here so the current verification scripts do not depend on the
development-history scripts archived in validation/archive/.
"""

import os
from concurrent.futures import ProcessPoolExecutor
from math import comb

# One BLAS/OpenMP thread per process: the classifier otherwise starts ~16
# OpenMP threads and the pool workers fight over the cores. Scripts import
# this module before numpy, so the setting takes hold; workers inherit it.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

# USDA classes in texture-triangle order, used for per-class tables.
ORDER = ["sand", "loamy sand", "sandy loam", "loam", "silt", "silt loam",
         "sandy clay loam", "clay loam", "silty clay loam", "sandy clay",
         "silty clay", "clay"]

# The reference columns the kNN needs.
COLS = ["layer_id", "profile_id", "texture_class", "alpha_kpa", "n", "thetar",
        "thetas", "sand", "silt", "clay", "ksat_cmh", "depth_cm"]

# The four aggregated texture groups.
GROUP = {
    "sand": "Sandy", "loamy sand": "Sandy",
    "sandy loam": "Loamy", "loam": "Loamy",
    "sandy clay loam": "Loamy", "clay loam": "Loamy",
    "silt": "Silty", "silt loam": "Silty", "silty clay loam": "Silty",
    "sandy clay": "Clayey", "silty clay": "Clayey", "clay": "Clayey",
}
GROUP_ORDER = ["Sandy", "Loamy", "Silty", "Clayey"]

# Bounding box used for "European" soils.
EU_BBOX = dict(lat=(34, 72), lon=(-25, 45))


def n_jobs():
    """Worker processes for map_folds: every core, unless SWCC_JOBS says."""
    return int(os.environ.get("SWCC_JOBS", os.cpu_count()))


def map_folds(fn, folds, *args, threads=1, jobs=None):
    """[fn(fold, *args) for fold in folds], one worker process per fold.

    Folds are independent, so running them side by side uses every core far
    better than threading one model at a time: a gradient-boosted model on a
    few features keeps 16 threads only partly busy. Each worker is limited to
    `threads` BLAS/OpenMP threads so the workers do not oversubscribe the
    CPU. Pass large arrays through `args` rather than a closure: joblib
    memory-maps big numpy arguments instead of copying them to every worker.
    Results come back in fold order. SWCC_JOBS=1 runs serially, for debugging.
    """
    from joblib import Parallel, delayed, parallel_config
    jobs = jobs or n_jobs()
    if jobs == 1:
        return [fn(f, *args) for f in folds]
    with parallel_config(backend="loky", inner_max_num_threads=threads):
        return Parallel(n_jobs=jobs)(delayed(fn)(f, *args) for f in folds)


DATA = {}   # in a pool worker: the tables pool() sent


def _init(data):
    DATA.update(data)


def pool(data):
    """A process pool on every core. `data`, a dict of the large tables,
    reaches each worker once, as DATA, instead of once per task."""
    return ProcessPoolExecutor(n_jobs(), initializer=_init, initargs=(data,))


def merge(out, part, sel):
    """Copy rows `sel` of a task's result into `out` (arrays, or dicts of
    them, nested as the script keeps them)."""
    if isinstance(out, dict):
        for k in out:
            merge(out[k], part[k], sel)
    else:
        out[sel] = part[sel]


def europe_mask(df):
    return (df.lat.between(*EU_BBOX["lat"]) & df.lon.between(*EU_BBOX["lon"])
            ).fillna(False)


def mcnemar(a_ok, b_ok):
    """Exact two-sided McNemar p-value for paired correct/incorrect arrays."""
    b = int(np.sum(a_ok & ~b_ok))
    c = int(np.sum(~a_ok & b_ok))
    n = b + c
    if n == 0:
        return 1.0
    lo = min(b, c)
    return min(1.0, sum(comb(n, i) for i in range(lo + 1)) / 2 ** n * 2)


def prf(truth, preds):
    """Macro recall / precision / F1 over the classes present in truth."""
    recs, precs, f1s = [], [], []
    for c in ORDER:
        m, pm = truth == c, preds == c
        if not m.any():
            continue
        tp = int((m & pm).sum())
        r = tp / m.sum()
        p = tp / pm.sum() if pm.sum() else 0.0
        recs.append(r); precs.append(p)
        f1s.append(2 * p * r / (p + r) if (p + r) else 0.0)
    return np.mean(recs), np.mean(precs), np.mean(f1s)
