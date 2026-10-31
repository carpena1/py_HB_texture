"""Donana sands through the Laikipia shrinkage conversion.

Same steps as the published Laikipia volumetric series (Willard et al. 2026,
HydroShare): w = theta / BD_sample, BD(w) = min(cap, 1 / (1/2.65 + w)),
theta' = w * BD(w).
BD_sample = 1.70 (mean of the 47 Donana rings, 1.59-1.80). Two caps: the
Laikipia/Tomejil 1.6 (a vertisol's shrinkage limit) and the sands' own
measured density 1.70. Panel c compares theta'/theta with the 70 published
Laikipia curves (cap 1.6, KALRO density; see test_caps.samples).
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
os.chdir(ROOT)
sys.argv = ["x"]
import prepare_andalusia as pa

BD, RHO_S = pa.DONANA_BD, 2.65
OUT = os.path.join(HERE, "donana_shrinkage_check")
GREY, BROWN, LIGHT = "#7f8c8d", "#7b4a35", "#c08a66"
KPA_PER_CM = 0.0980665


def convert(theta, cap):
    w = theta / BD
    bd = np.minimum(cap, 1.0 / (1.0 / RHO_S + w))
    return w * bd, bd, 1.0 / (1.0 / RHO_S + w)


def pf(h_kpa):
    with np.errstate(divide="ignore"):
        return np.where(h_kpa > 0, np.log10(h_kpa / KPA_PER_CM), -0.5)


curves = pa._donana_curves()
site, key = "AC", "020"
depth, h, th = curves[(site, key)]
o = np.argsort(h)
h, th = h[o], th[o]
x = pf(h)
own, bd_own, line = convert(th, BD)
lai, bd_lai, _ = convert(th, 1.6)

fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(16, 5.2))
fig.suptitle("Doñana sand (El Acebuche, 20 cm): the shrinking-density conversion that changed the "
             "Laikipia clays leaves these curves unchanged", fontsize=12.5)

# a: water content
a1.plot(x, th, "o-", color=GREY, lw=3, ms=8, label=f"as measured, one density ({BD:.2f})")
a1.plot(x, own, "o-", color=BROWN, lw=1.6, ms=4,
        label="shrinking density, capped at its own 1.70")
a1.plot(x, lai, "s--", color=LIGHT, lw=1.4, ms=4,
        label="shrinking density, Laikipia cap 1.6")
a1.axhline(1 - BD / RHO_S, color=GREY, ls=":", lw=1)
a1.text(5.2, 1 - BD / RHO_S + 0.006, f"porosity at 1.70: {1 - BD / RHO_S:.3f}", ha="right",
        fontsize=8.5, color=GREY)
a1.set_xlabel("pF  (log10 suction in cm; h = 0 drawn at −0.5)")
a1.set_ylabel("volumetric water content θ (m³/m³)")
a1.set_ylim(-0.01, 0.42)
a1.legend(fontsize=7.6, loc="center left", bbox_to_anchor=(0.40, 0.55), frameon=False, handlelength=1.8)
a1.set_title("a  retention curve", loc="left", fontsize=10.5)

# b: bulk density
a2.axhspan(1.59, 1.80, color=GREY, alpha=0.12, lw=0)
a2.text(-0.4, 1.795, "range of the 47 rings (1.59–1.80)", fontsize=8.3, color=GREY, va="top")
a2.plot(x, np.full_like(x, BD), "o-", color=GREY, lw=2, ms=6, label="one density (1.70)")
a2.plot(x, np.minimum(line, 2.0), ":", color=BROWN, lw=1.2,
        label="normal-shrinkage line 1/(1/2.65 + w), uncapped")
a2.plot(x, bd_own, "o-", color=BROWN, lw=1.6, ms=4, label="BD(w) capped at 1.70")
a2.plot(x, bd_lai, "s--", color=LIGHT, lw=1.4, ms=4, label="BD(w) capped at 1.6 (Laikipia rule)")
a2.set_ylim(1.0, 2.0)
a2.set_xlabel("pF")
a2.set_ylabel("bulk density (g/cm³)")
a2.legend(fontsize=8.3, loc="center right", bbox_to_anchor=(1.0, 0.36), frameon=False)
a2.set_title("b  density used for the conversion", loc="left", fontsize=10.5)
w_sat = (1 - BD / RHO_S) / BD
a2.text(-0.55, 1.03, f"saturated sand: w = {w_sat:.3f}, already on the shrinkage line\n"
        f"(Laikipia clays at saturation: w ≈ 0.5–0.7, density 0.9–1.2)",
        ha="left", fontsize=8.3, color="#444444")

# c: theta'/theta, all Donana layers against all Laikipia samples
sys.path.insert(0, HERE)
import test_caps as tc                                       # noqa: E402
lk = pd.concat([pd.DataFrame({"SampleID": sid, "pF": np.round(np.where(
    h_ > 0, np.log10(np.maximum(h_, 1e-9) / 0.0980665), 0.0), 2),
    "ratio": np.minimum(1.6, 1.0 / (1.0 / 2.65 + w_)) / bd1_})
    for sid, h_, w_, bd1_, _, _ in tc.samples("laikipia")], ignore_index=True)
q = lk.groupby("pF").ratio.quantile([0.25, 0.5, 0.75]).unstack()
a3.fill_between(q.index, q[0.25], q[0.75], color=BROWN, alpha=0.18, lw=0)
a3.plot(q.index, q[0.5], "o-", color=BROWN, lw=2, ms=5,
        label=f"Laikipia vertic clays ({lk.SampleID.nunique()} samples): median, 25–75 %")
first = True
r_own_all, r_lai_all = [], []
for (s_, k_), (_, hh, tt) in curves.items():
    ok = tt > 0.005
    hh, tt = hh[ok], tt[ok]
    ro = convert(tt, BD)[0] / tt
    rl = convert(tt, 1.6)[0] / tt
    r_own_all += list(ro); r_lai_all += list(rl)
    a3.plot(pf(hh), ro, "-", color=GREY, lw=0.8, alpha=0.7,
            label=f"Doñana sands ({len(curves)} layers), cap 1.70" if first else None)
    a3.plot(pf(hh), rl, "--", color=LIGHT, lw=0.8, alpha=0.7,
            label="Doñana sands, Laikipia cap 1.6" if first else None)
    first = False
a3.axhline(1.0, color="#999999", lw=0.8)
a3.set_ylim(0.85, 1.75)
a3.set_xlabel("pF")
a3.set_ylabel("θ converted / θ as measured")
a3.legend(fontsize=8.3, loc="upper left", frameon=False)
a3.set_title("c  how much the conversion changes θ", loc="left", fontsize=10.5)
a3.text(5.25, 1.10, f"Doñana, cap 1.70: {np.min(r_own_all):.3f}–{np.max(r_own_all):.3f}\n"
        f"Laikipia at pF 4.2: median {q.loc[4.2, 0.5]:.2f}", ha="right", fontsize=8.3, color="#444444")
for ax in (a1, a2, a3):
    ax.set_xlim(-0.7, 5.4)
    ax.grid(True, color="#e6e6e6", lw=0.6)
    ax.set_axisbelow(True)
fig.tight_layout(rect=(0, 0, 1, 0.94))
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}.{ext}", dpi=200, bbox_inches="tight")
print(f"wrote {OUT}.png/.pdf")
print(f"example layer: theta change with cap 1.70: {np.nanmax(np.abs(own / np.where(th > 0.005, th, np.nan) - 1)) * 100:.1f} % max;"
      f" with cap 1.6: {np.nanmedian(lai[th > 0.005] / th[th > 0.005]):.3f}")
print(f"all Donana points: cap 1.70 ratio {np.min(r_own_all):.4f}-{np.max(r_own_all):.4f}; cap 1.6 ratio "
      f"{np.min(r_lai_all):.4f}-{np.max(r_lai_all):.4f}")
print(q.round(3).to_string())
