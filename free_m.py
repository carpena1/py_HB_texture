"""Unconstrained-m van Genuchten fitting, shared by the prepare_*.py scripts.

The Mualem simplification ties the two shape parameters together as
m = 1 - 1/n. That is a convenience for deriving conductivity, not a property
of soils: n and m govern different parts of the retention curve, and the dry
limb -- the part most closely tied to texture -- is controlled mainly by m.
Freeing m therefore gives the texture inference a coordinate the constrained
fit cannot express.

Every reference table carries both parameter sets:

    thetar,   thetas,   alpha_kpa,   n         Mualem-constrained (m = 1-1/n)
    thetar_m, thetas_m, alpha_kpa_m, n_m,  m   unconstrained

Where a database publishes only fitted parameters and no measured points, the
unconstrained columns are copied from the constrained ones with m = 1 - 1/n,
which is exactly what those published parameters assume.
"""

import numpy as np

import swcc_texture as st

MIN_PTS_FREE = 6          # 5 parameters plus one degree of freedom
FREE_COLS = ["thetar_m", "thetas_m", "alpha_kpa_m", "n_m", "m", "rmse_m"]


def mualem_fallback(thetar, thetas, alpha_kpa, n, rmse=np.nan):
    """The unconstrained columns implied by a Mualem-constrained fit."""
    return dict(thetar_m=thetar, thetas_m=thetas, alpha_kpa_m=alpha_kpa,
                n_m=n, m=1.0 - 1.0 / n, rmse_m=rmse)


def fit_free_m(h, theta, fallback=None):
    """Fit the five-parameter vG form to measured (h [kPa], theta) points.

    Returns the FREE_COLS dict. Falls back to `fallback` -- normally the
    result of mualem_fallback() -- when there are too few points or the fit
    does not converge.
    """
    h = np.asarray(h, dtype=float)
    theta = np.asarray(theta, dtype=float)
    ok = np.isfinite(h) & np.isfinite(theta) & (h >= 0) & (theta > 0)
    h, theta = h[ok], theta[ok]
    if len(h) < MIN_PTS_FREE:
        return fallback
    try:
        popt, _ = st.fit_vg(h, theta, free_m=True)
    except Exception:
        return fallback
    tr, ts, la, ln1, m = popt
    a, n = 10.0 ** la, 1.0 + 10.0 ** ln1
    rmse = float(np.sqrt(np.mean((st.vg_theta(h, tr, ts, a, n, m) - theta) ** 2)))
    return dict(thetar_m=tr, thetas_m=ts, alpha_kpa_m=a, n_m=n, m=m,
                rmse_m=rmse)
