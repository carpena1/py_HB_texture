"""Two vertisols, two regions, one conversion: Laikipia (top) and Tomejil
(bottom), each converted from one density to a density that follows water
content (normal shrinkage, BD(w) = min(1.75, 1/(1/2.65 + w))); top right, how
much the conversion changes theta across all samples of both.

Laikipia: data/willard_raw/corrected/SWCC_shrinkage_summary.csv (one density
= each sample's own, behind the corrected porosity; pF 0 is saturation, drawn
at h = 0). Tomejil: prepare_andalusia._tomejil_curves (gravimetric; one
density = treatment mean, 1.085 CT / 1.125 DD).
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

ONE = {"CT": 1.085, "DD": 1.125}
CAP = 1.75     # median oven-dry clod density, KSSL Vertisols with clay 45-65 % and BD33 < 1.25
OUT = os.path.join(HERE, "tomejil_shrinkage_vs_laikipia")
GREY, BROWN, LAI = "#7f8c8d", "#7b4a35", "#c08a66"
LK_EX, TJ_EX = "3LT07414", "DD30"
XL = (-0.8, 6.7)


def pf_cm(h_cm):
    with np.errstate(divide="ignore"):
        return np.where(h_cm > 0, np.log10(h_cm), -0.5)


def lk_pf(p):
    return np.where(np.asarray(p) == 0, -0.5, p)          # pF 0 = saturation


lk = pd.read_csv("data/willard_raw/corrected/SWCC_shrinkage_summary.csv")
lk = lk[lk.w.notna()].copy()
# the summary's shrinkage columns were computed with a 1.6 cap: recompute with CAP
lk["bd_shrink"] = np.minimum(CAP, 1.0 / (1.0 / 2.65 + lk.w))
lk["shrinkage"] = 100 * lk.w * lk.bd_shrink
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
    a.set_ylim(0, 0.65)
    a.set_ylabel("volumetric water content θ (m³/m³)")
    a.legend(fontsize=8.3, loc="lower left", frameon=False)
    a.set_title(title, loc="left", fontsize=10.5)


def density_panel(a, x, bd1, bds, lab1, title):
    a.plot(x, np.full_like(x, bd1, dtype=float), "o-", color=GREY, lw=2, ms=5, label=lab1)
    a.plot(x, bds, "o-", color=BROWN, lw=2, ms=5, label=f"BD(w) = min({CAP:g}, 1/(1/2.65 + w))")
    a.axhline(CAP, color="#999999", ls="--", lw=1)
    a.text(XL[0] + 0.2, CAP - 0.02, f"cap {CAP:g} (shrinkage limit, w = {1 / CAP - 1 / 2.65:.3f})",
           fontsize=8.3, color="#666666", va="top")
    a.set_ylim(0.95, 1.85)
    a.set_ylabel("bulk density (g/cm³)")
    a.legend(fontsize=8.3, loc="center right", bbox_to_anchor=(1.0, 0.42), frameon=False)
    a.set_title(title, loc="left", fontsize=10.5)


# top row: Laikipia
xl = lk_pf(ex.pF.to_numpy())
curve_panel(ax[0, 0], xl, ex.corrected / 100, ex.shrinkage / 100,
            f"one density ({ex.bd_sample.iloc[0]:.3f} g/cm³, the sample's own)",
            f"a  Laikipia (Kenya), sample WL_{LK_EX}: retention curve")
density_panel(ax[0, 1], xl, ex.bd_sample.iloc[0], ex.bd_shrink.to_numpy(),
              f"one density ({ex.bd_sample.iloc[0]:.3f})", "b  Laikipia: density used for the conversion")
ax[0, 0].text(np.log10(15000) + 0.05, 0.63, "1500 kPa\n(the tool's\ndry limit)", fontsize=7.8,
              color="#888888", va="top")

# bottom row: Tomejil
xt = pf_cm(h)
bd1 = ONE[treat]
curve_panel(ax[1, 0], xt, w * bd1, w * pa.shrinkage_bd(w, CAP),
            f"one density ({bd1:.3f} g/cm³, treatment mean)",
            f"d  Tomejil (Andalusia), ring {TJ_EX}: retention curve")
density_panel(ax[1, 1], xt, bd1, pa.shrinkage_bd(w, CAP), f"one density ({bd1:.3f})",
              "e  Tomejil: density used for the conversion")

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
a.set_ylim(0.9, 1.7)
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
    f"   BD(w) = min({CAP:g}, 1 / (1/2.65 + w))   normal shrinkage:\n"
    "      the soil loses volume equal to the water it loses,\n"
    f"      until the shrinkage limit (BD {CAP:g}, w = {1 / CAP - 1 / 2.65:.3f};\n"
    "      median oven-dry clod density of KSSL Vertisols)\n"
    "   θ = w · BD(w)\n\n"
    "One density:\n"
    "   Laikipia — each sample's own, behind the corrected porosity\n"
    "   Tomejil — no ring density; the treatment mean,\n"
    "      1.085 (conventional tillage) or 1.125 (direct drill)\n"
    "      (Ordóñez Fernández et al. 2007); curves gravimetric\n\n"
    "Both soils start near 1.1–1.2 g/cm³ when saturated, so the\n"
    "conversion lifts the dry end by about half. Beyond pF ~3.8 most\n"
    "Tomejil rings reach the cap, the ratio then cap / treatment density\n"
    f"   ({CAP:g} / 1.125 = {CAP / 1.125:.2f};  {CAP:g} / 1.085 = {CAP / 1.085:.2f}).\n\n"
    "The agreement is partly built in (same formula and cap).\n"
    "Independent support: with the converted curves the tool\n"
    "reads both as clay — Tomejil 54 % clay predicted (56 %\n"
    "measured), Laikipia 52 % (52 %). A measured shrinkage\n"
    "curve or COLE would test the model directly.")
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
