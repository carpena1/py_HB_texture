"""Saturated conductivity from the retention curve alone, by capillary theory.

The soil is treated as a bundle of capillaries: the retention curve gives the
pore-size distribution through the capillary equation r = 2 gamma / (rho g h),
each pore conducts by Hagen-Poiseuille, and Ks is the sum over the pores
(Hillel 1980, ch. 8-A, 8-B, 8-I and 9-G). Unlike the kNN, this uses no
reference data at all, so it is an independent opinion on the same curve.

Two forms, both absolute (no measured Ks to scale them):

  marshall   Marshall (1958) as Hillel gives it: the water content range is cut
             into m equal increments, each increment's midpoint suction gives a
             pore radius, and
                 k = (eps^2 / 8 m^2) sum_i (2i - 1) r_i^2,   Ks = k rho g / eta
             with the radii in decreasing order, so the widest pore carries
             weight 1 and the narrowest 2m-1. This is the classical pairing
             statistic of Childs and Collis-George (1950).

  peters     Peters et al. (2023, HESS 27:1565), the same capillary bundle with
             the conductivity written so that no matching factor is needed:
                 Ks = beta tau_s (theta_s - theta_r)^2 alpha^2
             for the constrained van Genuchten curve (m = 1 - 1/n), with
             beta = 3.04e-4 m3/s at 20 C and a saturated tortuosity tau_s that
             they find is nearly soil-independent: 0.062 for this curve (0.1 as
             a general default). alpha is in m-1.

Both return cm/h. TAU_S and the physical constants are the published values.

Two choices were made on the reference's measured Ks (verify_ks_physical.py),
and the tool uses both:

  AIR_ENTRY_CM  marshall_ks counts no pore wider than the one that empties at
                this suction. Uncapped, a van Genuchten curve with n < 2 puts
                ever wider pores near saturation, the widest increment
                dominates the sum, and the result runs high by a factor of
                six. The cap was chosen out of fold from 2 to 200 cm, and
                every fold chose 50 cm or more (50 cm in 13 of 19). 50 cm is
                a pore about 60 um across, near the 50 um lower limit of
                Greenland's (1977) transmission pores.
  MATCHING_FACTOR  the classical scaling to a measured Ks (Childs and
                Collis-George 1950; Jackson 1972), since the curve fixes the
                shape of the pore-size distribution but not its absolute
                level. Applied by the caller, which reports the value with and
                without it.
"""

import numpy as np

# Water at 20 C, SI.
SURFACE_TENSION = 0.0728        # N/m
VISCOSITY = 1.002e-3            # N s/m2
DENSITY = 998.2                 # kg/m3
GRAVITY = 9.81                  # m/s2
BETA = 3.04e-4                  # m3/s, Peters et al. (2023) eq. 12
TAU_S = 0.062                   # their constrained-vG value; 0.1 in general
KPA_PER_M = 9.80665             # head of water: 1 m = 9.80665 kPa
M_S_TO_CM_H = 100.0 * 3600.0
N_INCREMENTS = 100              # Marshall's m
AIR_ENTRY_CM = 50.0             # see the module docstring

# Median over the laboratories of each one's median marshall / measured Ks
# (with the cap), so each source has one vote and the largest (Florida, half
# the layers) does not set it. 1.085 on the merged reference (23 sources with
# a measured Ks), 1.087 on the public one (13); the sources themselves range
# from 0.15 to 200. Without the cap it would be 3.9, the sources ranging from
# 0.8 to 8,000 (verify_ks_physical.py).
MATCHING_FACTOR = 1.09


def _to_cm_h(k_m_s):
    return np.asarray(k_m_s, float) * M_S_TO_CM_H


def peters_ks(thetar, thetas, alpha_kpa, tau_s=TAU_S):
    """Ks (cm/h) from Peters et al. (2023), constrained van Genuchten."""
    alpha_m = np.asarray(alpha_kpa, float) * KPA_PER_M          # 1/kPa -> 1/m
    ks = BETA * tau_s * (np.asarray(thetas, float)
                         - np.asarray(thetar, float)) ** 2 * alpha_m ** 2
    return _to_cm_h(ks)


def marshall_ks(thetar, thetas, alpha_kpa, n, m=None, m_inc=N_INCREMENTS,
                h_min_cm=AIR_ENTRY_CM):
    """Ks (cm/h) from Marshall (1958) on the fitted van Genuchten curve.

    The pore radii come from the suction at the midpoint of each of m_inc
    equal water-content increments between thetar and thetas. m defaults to
    the constrained 1 - 1/n; pass it when the curve was fitted with m free.

    h_min_cm caps the pore radius at that suction's, an air-entry value: for
    n < 2 the van Genuchten curve puts ever wider pores near saturation, and
    the widest increment can dominate the sum (Vogel et al. 2001; Ippisch et
    al. 2006). None removes the cap.
    """
    thetar, thetas = np.atleast_1d(thetar).astype(float), np.atleast_1d(thetas).astype(float)
    alpha_m = np.atleast_1d(alpha_kpa).astype(float) * KPA_PER_M
    n = np.atleast_1d(n).astype(float)
    mm = 1.0 - 1.0 / n if m is None else np.atleast_1d(m).astype(float)

    # Midpoint effective saturation of each increment, widest pore first.
    i = np.arange(1, m_inc + 1)
    se = ((m_inc - i + 0.5) / m_inc)[None, :]
    # Inverse van Genuchten: h (m of water) at that saturation.
    h_m = (se ** (-1.0 / mm[:, None]) - 1.0) ** (1.0 / n[:, None]) / alpha_m[:, None]
    if h_min_cm is not None:
        h_m = np.maximum(h_m, h_min_cm / 100.0)
    r = 2.0 * SURFACE_TENSION / (DENSITY * GRAVITY * h_m)       # capillary, m

    eps = (thetas - thetar)[:, None]        # the pore space the curve drains
    k = (eps ** 2 / (8.0 * m_inc ** 2)) * ((2 * i - 1)[None, :] * r ** 2)
    ks = k.sum(axis=1) * DENSITY * GRAVITY / VISCOSITY
    return _to_cm_h(ks)
