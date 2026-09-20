"""Improved van Genuchten (IvG) retention model, for testing as an
alternative to the vG parameters the tool matches on.

Ghorbani, A., Babaeian, E., Sadeghi, M., Durner, W., Jones, S.B. and van
Genuchten, M.Th. (2025). An improved van Genuchten soil water characteristic
model to account for surface adsorptive forces. Journal of Hydrology 661,
133692, eq. (6):

    theta = thetas Se + thetar ln(h0/h) (1 - Se),   Se = [1 + (alpha h)^n]^-m

with h0 = 10^6.8 cm, where soils are oven dry and theta = 0. The second term
replaces vG's constant residual water content with a log-linear dry limb, so
thetar is a fitting parameter rather than a residual water content.

Because the model forces theta = 0 at h0 by construction, adding an
(h0, 0) point to the data changes nothing; what the dry end does depend on is
the bound on thetar. Two are offered, and verify_ivg.py compares them:

    bound="thetas20"   thetar <= thetas/20, the paper's constraint (8), which
                       also keeps the curve from rising above thetas on its
                       way to zero at h0
    bound="abs"        thetar <= 0.03, as in the Arizona IvG fits
                       (IvG_curves/AZ_SoilData_Basic_SWC_vG_IvG_parameters.xlsx)

h is in kPa, as in swcc_texture.
"""

import numpy as np
from scipy.optimize import curve_fit

import swcc_texture as st

H0_KPA = 10.0 ** 6.8 * 0.0980665      # 10^6.8 cm, about 6.2e5 kPa
THETAR_MAX = 0.03                      # bound="abs"
C_MAX = 20.0                           # bound="thetas20", eq. (8)


def ivg_theta(h, thetar, thetas, alpha, n, m=None):
    """IvG retention at h (kPa >= 0); m defaults to 1 - 1/n."""
    if m is None:
        m = 1.0 - 1.0 / n
    h = np.asarray(h, dtype=float)
    se = (1.0 + (alpha * h) ** n) ** (-m)
    hp = np.where(h > 0, h, 1.0)
    log_term = np.where(h > 0, np.log(H0_KPA / np.minimum(hp, H0_KPA)), 0.0)
    return thetas * se + thetar * log_term * (1.0 - se)


def _transformed(h, thetar, thetas, la, ln1, m=None):
    return ivg_theta(h, thetar, thetas, 10.0 ** la, 1.0 + 10.0 ** ln1, m)


def _transformed_frac(h, r, thetas, la, ln1, m=None):
    # thetar = r thetas, so the bound thetar <= thetas / C_MAX is a box.
    return _transformed(h, r * thetas, thetas, la, ln1, m)


def fit_ivg(h, theta, free_m=False, bound="abs"):
    """Fit IvG to (h [kPa], theta) points.

    Returns (thetar, thetas, alpha_kpa, n, m, rmse), with m = 1 - 1/n unless
    free_m, and thetar bounded as `bound` says (see the module docstring). Bounds follow swcc_texture.fit_vg. The log term gives the fit
    several local minima (thetar at its bound with a low n, or near zero with
    a steep n), so it starts from a grid of thetar and n and keeps the best.
    """
    h = np.asarray(h, dtype=float)
    theta = np.asarray(theta, dtype=float)
    npar = 5 if free_m else 4
    if len(h) < npar + 1:
        raise ValueError(f"need at least {npar + 1} points")
    tmax = theta.max()
    hpos = h[h > 0]
    frac = bound == "thetas20"
    top = 1.0 / C_MAX if frac else THETAR_MAX
    lb = [0.0, 0.5 * tmax, -4.0, np.log10(0.01)]
    ub = [top, 1.0, 1.5, np.log10(10.0)]
    g = _transformed_frac if frac else _transformed
    if free_m:
        f = g
        lb, ub = lb + [st.M_BOUNDS[0]], ub + [st.M_BOUNDS[1]]
    else:
        def f(h, tr, ts, la, ln1):
            return g(h, tr, ts, la, ln1)
    best = None
    for tr0 in (0.1, 0.5, 0.9):
        for ln10 in (np.log10(0.3), np.log10(1.0), np.log10(3.0)):
            p0 = [tr0 * top, tmax, np.log10(1.0 / np.median(hpos)), ln10]
            if free_m:
                p0 = p0 + [1.0 - 1.0 / (1.0 + 10.0 ** ln10)]
            try:
                popt, _ = curve_fit(f, h, theta, p0=np.clip(p0, lb, ub),
                                    bounds=(lb, ub), maxfev=20000)
            except RuntimeError:
                continue
            sse = float(np.sum((f(h, *popt) - theta) ** 2))
            if best is None or sse < best[0]:
                best = (sse, popt)
    if best is None:
        raise RuntimeError("IvG fit did not converge from any start")
    popt = best[1]
    tr, ts, la, ln1 = popt[:4]
    if frac:
        tr = tr * ts
    alpha, n = 10.0 ** la, 1.0 + 10.0 ** ln1
    m = popt[4] if free_m else 1.0 - 1.0 / n
    rmse = float(np.sqrt(np.mean((ivg_theta(h, tr, ts, alpha, n, m)
                                  - theta) ** 2)))
    return tr, ts, alpha, n, m, rmse
