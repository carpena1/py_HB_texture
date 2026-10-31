"""Two vertisols, two regions, one conversion: Laikipia (top) and Tomejil
(bottom), each converted from one density to a density that follows water
content (normal shrinkage, BD(w) = min(cap, 1/(1/2.65 + w))); top right, how
much the conversion changes theta across all samples of both.

Laikipia: the published gravimetric curves (Willard et al. 2026, HydroShare,
doi:10.4211/hs.99080b48206e4b2c82e716516d463a61; 70 samples, the dispersed
ones left out), one density = the KALRO density of the core, cap 1.6 as in the
published volumetric series; pF 0 is saturation, drawn at h = 0. Tomejil:
prepare_andalusia._tomejil_curves (gravimetric; one density = treatment mean,
1.085 CT / 1.125 DD), cap 1.75.
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
import test_caps as tc

ONE = {"CT": 1.085, "DD": 1.125}
CAP_LK = 1.6   # Laikipia: the cap of the published volumetric series
CAP = 1.75     # Tomejil: median oven-dry clod density, KSSL Vertisols with clay 45-65 % and BD33 < 1.25
OUT = os.path.join(HERE, "tomejil_shrinkage_vs_laikipia")
GREY, BROWN, LAI = "#7f8c8d", "#7b4a35", "#c08a66"
LK_EX, TJ_EX = "3LT07414", "DD30"
XL = (-0.8, 6.7)


def pf_cm(h_cm):
    with np.errstate(divide="ignore"):
        return np.where(h_cm > 0, np.log10(h_cm), -0.5)


def lk_pf(p):
    return np.where(np.asarray(p) == 0, -0.5, p)          # pF 0 = saturation


rows = []
for sid, h_k, w_k, bd1_k, _, _ in tc.samples("laikipia"):
    pf_k = np.where(h_k > 0, np.log10(np.maximum(h_k, 1e-9) / 0.0980665), 0.0)
    bds_k = np.minimum(CAP_LK, 1.0 / (1.0 / 2.65 + w_k))
    rows.append(pd.DataFrame({"SampleID": sid, "pF": np.round(pf_k, 2), "w": w_k,
                              "bd_sample": bd1_k, "bd_shrink": bds_k}))
lk = pd.concat(rows, ignore_index=True)
lk["one"] = lk.w * lk.bd_sample
lk["shrinkage"] = lk.w * lk.bd_shrink
lk["ratio"] = lk.bd_shrink / lk.bd_sample
ex = lk[lk.SampleID == LK_EX].sort_values("pF")
tj = pa._tomejil_curves()
treat, _, h, w, _ = tj[TJ_EX]
o = np.argsort(h)
h, w = h[o], w[o]

fig, ax = plt.subplots(2, 3, figsize=(16, 9.6))
fig.suptitle("Two vertisols, two regions: converting retention curves with a shrinking bulk density",
             fontsize=13.5, y=0.995)


def curve_panel(a, x, th1, ths, lab1, title, dry_limit=True):
    a.plot(x, th1, "o-", color=GREY, lw=2, ms=5, label=lab1)
    a.plot(x, ths, "o-", color=BROWN, lw=2, ms=5, label="density follows water content (normal shrinkage)")
    if dry_limit:
        a.axvline(np.log10(15000), color="#bbbbbb", lw=0.8, ls="--")
    a.set_ylim(0, 0.75)
    a.set_ylabel("volumetric water content θ (m³/m³)")
    a.legend(fontsize=8.3, loc="lower left", frameon=False)
    a.set_title(title, loc="left", fontsize=10.5)


def density_panel(a, x, bd1, bds, lab1, title, cap):
    a.plot(x, np.full_like(x, bd1, dtype=float), "o-", color=GREY, lw=2, ms=5, label=lab1)
    a.plot(x, bds, "o-", color=BROWN, lw=2, ms=5, label=f"BD(w) = min({cap:g}, 1/(1/2.65 + w))")
    a.axhline(cap, color="#999999", ls="--", lw=1)
    a.text(XL[0] + 0.2, cap - 0.02, f"cap {cap:g} (shrinkage limit, w = {1 / cap - 1 / 2.65:.3f})",
           fontsize=8.3, color="#666666", va="top")
    a.set_ylim(0.85, 1.85)
    a.set_ylabel("bulk density (g/cm³)")
    a.legend(fontsize=8.3, loc="center right", bbox_to_anchor=(1.0, 0.42), frameon=False)
    a.set_title(title, loc="left", fontsize=10.5)


# top row: Laikipia
xl = lk_pf(ex.pF.to_numpy())
curve_panel(ax[0, 0], xl, ex.one, ex.shrinkage,
            f"one density ({ex.bd_sample.iloc[0]:.3f} g/cm³, the core's KALRO density)",
            f"a  Laikipia (Kenya), sample WL_{LK_EX}: retention curve")
density_panel(ax[0, 1], xl, ex.bd_sample.iloc[0], ex.bd_shrink.to_numpy(),
              f"one density ({ex.bd_sample.iloc[0]:.3f})", "b  Laikipia: density used for the conversion",
              CAP_LK)
ax[0, 0].text(np.log10(15000) + 0.05, 0.73, "1500 kPa\n(the tool's\ndry limit)", fontsize=7.8,
              color="#888888", va="top")

# bottom row: Tomejil
xt = pf_cm(h)
bd1 = ONE[treat]
curve_panel(ax[1, 0], xt, w * bd1, w * pa.shrinkage_bd(w, CAP),
            f"one density ({bd1:.3f} g/cm³, treatment mean)",
            f"d  Tomejil (Andalusia), ring {TJ_EX}: retention curve")
density_panel(ax[1, 1], xt, bd1, pa.shrinkage_bd(w, CAP), f"one density ({bd1:.3f})",
              "e  Tomejil: density used for the conversion", CAP)

# top right: combined
a = ax[0, 2]
q = lk.groupby("pF").ratio.quantile([0.25, 0.5, 0.75]).unstack()
qx = lk_pf(q.index.to_numpy())
rows = []
for name, (t_, _, hh, ww, _) in tj.items():
    r = pa.shrinkage_bd(ww, CAP) / ONE[t_]
    s = np.argsort(hh)
    a.plot(pf_cm(hh)[s], r[s], "-", color=BROWN, lw=0.5, alpha=0.22)
    rows += list(zip(pf_cm(hh), r))
tr = pd.DataFrame(rows, columns=["pF", "r"])
tr["bin"] = (tr.pF / 0.5).round() * 0.5
tq = tr.groupby("bin").r.quantile([0.25, 0.5, 0.75]).unstack()
a.fill_between(tq.index, tq[0.25], tq[0.75], color=BROWN, alpha=0.18, lw=0)
a.plot(tq.index, tq[0.5], "o-", color=BROWN, lw=2, ms=4, label=f"Tomejil ({len(tj)} rings): median, 25–75 %")
a.fill_between(qx, q[0.25], q[0.75], color=LAI, alpha=0.25, lw=0)
a.plot(qx, q[0.5], "s-", color=LAI, lw=2, ms=4, label=f"Laikipia ({lk.SampleID.nunique()} samples): median, 25–75 %")
a.axhline(1.0, color="#999999", lw=0.8)
a.axvline(np.log10(15000), color="#bbbbbb", lw=0.8, ls="--")
a.set_ylim(0.8, 1.7)
a.set_ylabel("θ shrinking density / θ one density")
a.legend(fontsize=8.3, loc="upper left", frameon=False)
a.set_title("c  both soils: how much the conversion changes θ", loc="left", fontsize=10.5)
t25, t42 = np.interp(2.5, tq.index, tq[0.5]), np.interp(4.2, tq.index, tq[0.5])
a.text(XL[1] - 0.1, 0.93, f"median at pF 2.5: Tomejil {t25:.2f}, Laikipia {q.loc[2.5, 0.5]:.2f}\n"
       f"median at pF 4.2: Tomejil {t42:.2f}, Laikipia {q.loc[4.2, 0.5]:.2f}",
       ha="right", fontsize=8.3, color="#444444")

# bottom right: method
m = ax[1, 2]
m.axis("off")
m.set_title("f  the conversion", loc="left", fontsize=10.5)
txt = (
    "At every suction:\n"
    "   w = θ / BD (one density)         gravimetric water\n"
    "   BD(w) = min(cap, 1 / (1/2.65 + w))   normal shrinkage:\n"
    "      the soil loses volume equal to the water it loses,\n"
    "      until the shrinkage limit (the cap)\n"
    "   θ = w · BD(w)\n\n"
    "Cap and one density:\n"
    f"   Laikipia — cap {CAP_LK:g}, as in the published series (HydroShare);\n"
    "      the core's KALRO density; curves published gravimetric\n"
    f"   Tomejil — cap {CAP:g}, median oven-dry clod density of KSSL\n"
    "      Vertisols; the treatment mean, 1.085 (conventional tillage)\n"
    "      or 1.125 (direct drill; Ordóñez Fernández et al. 2007)\n\n"
    "Where the soil is wetter than the core was dense (swollen,\n"
    "BD(w) below the one density), the conversion lowers θ; as it\n"
    "dries and shrinks, it raises θ, by about half at the dry end.\n"
    "Beyond pF ~3.8 most rings and cores reach their cap.\n\n"
    "The agreement is partly built in (same formula).\n"
    "Independent support: with the converted curves the tool\n"
    "reads both soils as clay. A measured shrinkage curve or\n"
    "COLE would test the model directly.")
m.text(0.0, 0.97, txt, va="top", ha="left", fontsize=9, family="DejaVu Sans", color="#333333",
       transform=m.transAxes, linespacing=1.35)

for a in (ax[0, 0], ax[0, 1], ax[0, 2], ax[1, 0], ax[1, 1]):
    a.set_xlim(*XL)
    a.grid(True, color="#e6e6e6", lw=0.6)
    a.set_axisbelow(True)
    a.set_xlabel("pF  (log10 suction in cm; h = 0 drawn at −0.5)")
fig.tight_layout(rect=(0, 0, 1, 0.97))
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}.{ext}", dpi=200, bbox_inches="tight")
print(f"wrote {OUT}.png/.pdf")
