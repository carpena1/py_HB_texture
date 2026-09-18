"""Validation figures, drawn from figures/cache/validation_predictions.csv.

Called by make_validation_figures.py. Needs matplotlib, which the tool itself
does not. All figures are aggregates (rates, a confusion matrix, hexbin
densities), so no individual restricted sample is shown.
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                              # noqa: E402
from matplotlib.colors import LinearSegmentedColormap        # noqa: E402
from matplotlib.patches import Rectangle                      # noqa: E402
from scipy.stats import spearmanr                             # noqa: E402

from verify_ceiling import bayes_bracket                      # noqa: E402

OUT = os.path.join("figures", "validation")

# Light-mode chart palette (dataviz reference instance).
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
S1, S2 = "#2a78d6", "#eb6834"               # categorical slots 1 and 2
SEQ = LinearSegmentedColormap.from_list(
    "seq", [SURFACE, "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf",
            "#184f95", "#0d366b"])

# Classes ordered by texture group, so the groups form contiguous blocks.
GROUPED = ["sand", "loamy sand",
           "sandy loam", "loam", "sandy clay loam", "clay loam",
           "silt", "silt loam", "silty clay loam",
           "sandy clay", "silty clay", "clay"]

# verify_holdout.py: the 1,711 targets with the curve alone, and the tool as
# run with --depth and --bulk-density on the HOLDOUT_N targets that carry both
# (name, exact, group, exact with both, group with both).
HOLDOUT_N = 1659
HOLDOUT = [("Soil alone\n(the standard)", 41.2, 64.8, 45.0, 66.7),
           ("Whole profile", 39.7, 63.8, 41.4, 65.0),
           ("Whole data source", 30.5, 58.3, 31.9, 57.9)]


def style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 10, "text.color": INK, "axes.labelcolor": INK2,
        "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
        "axes.facecolor": SURFACE, "figure.facecolor": SURFACE,
        "savefig.facecolor": SURFACE, "axes.spines.top": False,
        "axes.spines.right": False, "xtick.color": MUTED, "ytick.color": MUTED,
        "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
        "legend.frameon": False})


def title(fig, head, sub):
    # Offsets in inches, so the spacing holds whatever the figure height.
    h = fig.get_figheight()
    fig.text(0.01, 1 - 0.10 / h, head, ha="left", va="top", fontsize=12.5,
             fontweight="bold", color=INK)
    fig.text(0.01, 1 - 0.40 / h, sub, ha="left", va="top", fontsize=9.5,
             color=INK2)


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=200,
                    bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print(f"  wrote {OUT}/{name}.png/.svg")


def hgrid(ax, ticks):
    ax.set_xticks(ticks)
    for t in ticks:
        ax.axvline(t, color=GRID, lw=0.8, zorder=0)
    ax.tick_params(axis="y", length=0)


# ---------------------------------------------------------------------------

def fig_limits(d, group):
    """Tool accuracy against chance and the Cover & Hart ceiling."""
    truth = d.texture_class.to_numpy()
    rows = []
    for label, n_cls, f in (("Exact USDA class", 12, lambda x: x),
                            ("Texture group", 4, lambda x: group[x])):
        t = np.array([f(x) for x in truth])
        acc = lambda col: np.mean(np.array([f(x) for x in d[col]]) == t) * 100
        nn = acc("nn_class")
        lo, hi = bayes_bracket(1 - nn / 100, n_cls)
        rows.append(dict(label=label, chance=100 / n_cls, nn=nn, hi=hi * 100,
                         base=acc("pred_base"), cov=acc("pred_cov")))

    fig, ax = plt.subplots(figsize=(8.6, 3.5))
    fig.subplots_adjust(top=0.72, bottom=0.2, left=0.2, right=0.97)
    for y, r in zip((1, 0), rows):
        ax.add_patch(Rectangle((r["nn"], y - 0.2), r["hi"] - r["nn"], 0.4,
                               color=GRID, zorder=1, lw=0))
        ax.plot([r["hi"]] * 2, [y - 0.26, y + 0.26], color=INK, lw=1.4, zorder=2)
        ax.text(r["hi"] + 0.8, y + 0.2, f"ceiling {r['hi']:.0f} %",
                fontsize=9, color=INK, va="center")
        ax.plot([r["chance"]] * 2, [y - 0.26, y + 0.26], color=MUTED, lw=1.4)
        ax.text(r["chance"], y - 0.36, f"chance {r['chance']:.0f} %",
                fontsize=8.5, color=INK2, ha="center", va="top")
        for key, c, dy in (("base", S1, -1), ("cov", S2, 1)):
            ax.plot(r[key], y, "o", ms=10, color=c, mec=SURFACE, mew=2, zorder=4)
            ax.text(r[key], y + 0.3 * dy, f"{r[key]:.0f} %", fontsize=9,
                    color=INK, ha="center", va="bottom" if dy > 0 else "top")
    ax.set_yticks([1, 0], [r["label"] for r in rows], fontsize=10)
    ax.set_xlim(0, 100); ax.set_ylim(-0.6, 1.6)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0f} %")
    hgrid(ax, [0, 25, 50, 75, 100])
    ax.spines["left"].set_visible(False)
    h = [plt.Line2D([], [], marker="o", ls="", ms=9, color=S1, mec=SURFACE, mew=2),
         plt.Line2D([], [], marker="o", ls="", ms=9, color=S2, mec=SURFACE, mew=2),
         Rectangle((0, 0), 1, 1, color=GRID)]
    ax.legend(h, ["retention curve only", "+ depth and bulk density",
                  "achievable range (nearest neighbour to Cover & Hart ceiling)"],
              loc="upper center", bbox_to_anchor=(0.42, -0.14), ncol=3,
              fontsize=8.5, handletextpad=0.4, columnspacing=1.2)
    title(fig, "How close the tool gets to the best achievable accuracy",
          f"Leave-one-soil-out, {len(d):,} soils (equal numbers per class). "
          "The grey band runs from plain nearest-neighbour matching to the\n"
          "Cover & Hart (1967) bound for any method using the four curve "
          "parameters.")
    save(fig, "v1_accuracy_vs_limits")
    return rows


def fig_confusion(d, group):
    """Where the errors go: row-normalised confusion matrix."""
    idx = {c: i for i, c in enumerate(GROUPED)}
    m = np.zeros((12, 12))
    for t, p in zip(d.texture_class, d.pred_cov):
        m[idx[t], idx[p]] += 1
    pct = m / m.sum(axis=1, keepdims=True) * 100

    fig, ax = plt.subplots(figsize=(8.2, 7.6))
    fig.subplots_adjust(top=0.86, bottom=0.2, left=0.2, right=0.93)
    im = ax.imshow(pct, cmap=SEQ, vmin=0, vmax=80)
    for i in range(12):
        for j in range(12):
            v = pct[i, j]
            if i == j or v >= 15:
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=8.5,
                        color="white" if v >= 45 else INK,
                        fontweight="bold" if i == j else "normal")
    for b in (1.5, 5.5, 8.5):                     # group block boundaries
        ax.axhline(b, color=INK2, lw=0.8); ax.axvline(b, color=INK2, lw=0.8)
    ax.set_xticks(range(12), GROUPED, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(12), GROUPED, fontsize=9)
    ax.set_xlabel("predicted class", color=INK2)
    ax.set_ylabel("true class", color=INK2)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("% of the true class", color=INK2)
    cb.outline.set_visible(False)
    cb.ax.tick_params(color=MUTED, labelcolor=INK2)
    title(fig, "Where the errors go",
          "Share of each true class assigned to each predicted class, "
          "leave-one-soil-out with depth and bulk density.\nDiagonal = correct. "
          "Lines separate the four texture groups; cells under 15 % are unlabelled.")
    save(fig, "v2_confusion_matrix")


def fig_recall(d):
    """Per-class recall, curve only vs with depth and bulk density."""
    rec = {c: (np.mean(d.pred_base[d.texture_class == c] == c) * 100,
               np.mean(d.pred_cov[d.texture_class == c] == c) * 100,
               int((d.texture_class == c).sum())) for c in GROUPED}
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    fig.subplots_adjust(top=0.85, bottom=0.14, left=0.2, right=0.9)
    ys = np.arange(12)[::-1]
    for y, c in zip(ys, GROUPED):
        a, b, n = rec[c]
        ax.plot([a, b], [y, y], color=AXIS, lw=2, zorder=1)
        ax.plot(a, y, "o", ms=9, color=S1, mec=SURFACE, mew=2, zorder=3)
        ax.plot(b, y, "o", ms=9, color=S2, mec=SURFACE, mew=2, zorder=3)
        ax.text(101, y, f"n={n}", fontsize=8, color=MUTED, va="center")
    ax.axvline(100 / 12, color=MUTED, lw=1.2)
    ax.text(100 / 12 + 1, 11.55, "chance", fontsize=8.5, color=INK2)
    for b in (9.5, 5.5, 2.5):
        ax.axhline(b, color=GRID, lw=0.8)
    ax.set_yticks(ys, GROUPED, fontsize=9.5)
    ax.set_xlim(0, 100); ax.set_ylim(-0.7, 11.9)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0f} %")
    hgrid(ax, [0, 25, 50, 75, 100])
    ax.spines["left"].set_visible(False)
    h = [plt.Line2D([], [], marker="o", ls="", ms=9, color=c, mec=SURFACE, mew=2)
         for c in (S1, S2)]
    ax.legend(h, ["retention curve only", "+ depth and bulk density"],
              loc="upper center", bbox_to_anchor=(0.45, -0.07), ncol=2, fontsize=9)
    title(fig, "Correct class, by true class",
          "Share of each class predicted exactly, leave-one-soil-out. "
          "Horizontal lines separate the texture groups.")
    save(fig, "v3_recall_by_class")


def ks_stats(obs, med):
    e = np.abs(np.log10(med) - np.log10(obs))
    return (10 ** np.median(e), np.mean(e <= 1) * 100,
            spearmanr(obs, med).statistic, len(e))


def fig_ks_scatter(d):
    """Predicted vs measured Ks as a density (nearly all soils with a
    measured Ks are undisturbed, so one panel)."""
    ok = d.ksat_cmh.gt(0) & d.ks_med.gt(0)
    sets = [("All soils with a measured Ks", d[ok])]
    fig, axes = plt.subplots(1, 1, figsize=(6.4, 5.6), squeeze=False)
    axes = axes[0]
    fig.subplots_adjust(top=0.8, bottom=0.12, left=0.14, right=0.84)
    lim = (-3.5, 2.5)
    xx = np.array(lim)
    for ax, (name, s) in zip(axes, sets):
        x, y = np.log10(s.ksat_cmh), np.log10(s.ks_med)
        ax.fill_between(xx, xx - 1, xx + 1, color=GRID, alpha=0.55, lw=0, zorder=0)
        ax.fill_between(xx, xx - np.log10(2), xx + np.log10(2), color=GRID,
                        lw=0, zorder=0)
        ax.plot(xx, xx, color=INK2, lw=1, zorder=1)
        hb = ax.hexbin(x, y, gridsize=26, extent=lim + lim, mincnt=1, cmap=SEQ,
                       vmin=0, zorder=2, linewidths=0.2, edgecolors=SURFACE)
        f, w10, rho, n = ks_stats(s.ksat_cmh.to_numpy(), s.ks_med.to_numpy())
        ax.text(0.97, 0.03, f"{name}\nn = {n}\ntypical error ×{f:.1f}\n"
                f"within 10×: {w10:.0f} %\nrank correlation ρ = {rho:.2f}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
                color=INK, linespacing=1.4, zorder=5,
                bbox=dict(facecolor=SURFACE, edgecolor="none", alpha=0.9, pad=4))
        ax.set_xlim(lim); ax.set_ylim(lim)
        ticks = [-3, -2, -1, 0, 1, 2]
        ax.set_xticks(ticks, [f"{10.0 ** t:g}" for t in ticks])
        ax.set_yticks(ticks, [f"{10.0 ** t:g}" for t in ticks])
        ax.set_xlabel("measured Ks (cm/h)")
        ax.set_aspect("equal")
    axes[0].set_ylabel("predicted Ks (cm/h)")
    cax = fig.add_axes([0.87, 0.22, 0.02, 0.45])
    cb = fig.colorbar(hb, cax=cax)
    cb.set_label("soils per cell", color=INK2)
    cb.outline.set_visible(False)
    cb.ax.tick_params(color=MUTED, labelcolor=INK2)
    title(fig, "Predicted against measured saturated conductivity",
          "Leave-one-soil-out. Line = perfect prediction; darker band =\n"
          "within a factor of 2, lighter band = within a factor of 10.")
    save(fig, "v4_ks_predicted_vs_measured")


def fig_ks_cdf(d):
    """Share of soils whose Ks is within a factor F of the measurement."""
    ok = d.ksat_cmh.gt(0) & d.ks_med.gt(0) & d.ks_cov_med.gt(0)
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    fig.subplots_adjust(top=0.8, bottom=0.14, left=0.1, right=0.8)
    F = np.logspace(0, 3, 300)
    # The two curves nearly coincide, so values go in one box and identity in
    # the legend rather than as labels that would overprint each other.
    rows = []
    for name, col, c in (("retention curve only", "ks_med", S1),
                         ("+ bulk density", "ks_cov_med", S2)):
        s = d[ok]
        e = np.abs(np.log10(s[col]) - np.log10(s.ksat_cmh)).to_numpy()
        share = [np.mean(e <= np.log10(f)) * 100 for f in F]
        ax.plot(F, share, color=c, lw=2, solid_capstyle="round", label=name)
        vals = [np.mean(e <= np.log10(f)) * 100 for f in (2, 10)]
        for f, v in zip((2, 10), vals):
            ax.plot(f, v, "o", ms=8, color=c, mec=SURFACE, mew=2, zorder=3)
        rows.append(f"{name:<22s}{vals[0]:4.0f} %{vals[1]:7.0f} %")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.22), fontsize=9)
    ax.text(0.97, 0.04, f"{'':<22s}{'2×':>6s}{'10×':>9s}\n" + "\n".join(rows),
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5,
            family="monospace", color=INK, linespacing=1.5)
    ax.set_xscale("log")
    ax.set_xlim(1, 1000); ax.set_ylim(0, 100)
    ax.set_xticks([1, 2, 5, 10, 100, 1000], ["1", "2", "5", "10", "100", "1000"])
    ax.yaxis.set_major_formatter(lambda v, _: f"{v:.0f} %")
    for t in (25, 50, 75, 100):
        ax.axhline(t, color=GRID, lw=0.8, zorder=0)
    ax.set_xlabel("within a factor of … of the measured Ks")
    ax.set_ylabel("share of soils")
    title(fig, "How close the Ks estimate is",
          "Share of soils whose predicted Ks lies within a given factor of "
          "the measurement,\nleave-one-soil-out, with and without bulk "
          "density in the neighbour search.")
    save(fig, "v5_ks_within_factor")


def fig_holdout():
    """Accuracy at each hold-out level (what is hidden with the test soil),
    with the curve alone and with depth and bulk density supplied."""
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    fig.subplots_adjust(top=0.78, bottom=0.2, left=0.22, right=0.95)
    ys = np.arange(len(HOLDOUT))[::-1] * 1.25
    h = 0.24
    bars = ((S1, 1.0), (S1, 0.45), (S2, 1.0), (S2, 0.45))
    for y, (name, *vals) in zip(ys, HOLDOUT):
        for i, (v, (c, a)) in enumerate(zip((vals[0], vals[2], vals[1], vals[3]), bars)):
            dy = (1.5 - i) * (h + 0.02)
            ax.barh(y + dy, v, height=h, color=c, alpha=a, zorder=2)
            ax.text(v + 1, y + dy, f"{v:.1f} %", va="center", fontsize=8.5, color=INK)
    ax.set_yticks(ys, [n for n, *_ in HOLDOUT], fontsize=9.5)
    ax.set_xlim(0, 100)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0f} %")
    hgrid(ax, [0, 25, 50, 75, 100])
    ax.spines["left"].set_visible(False)
    hd = [Rectangle((0, 0), 1, 1, color=c, alpha=a) for c, a in bars]
    ax.legend(hd, ["exact class, curve alone", "exact class, + depth & bulk density",
                   "texture group, curve alone", "texture group, + depth & bulk density"],
              loc="upper center", bbox_to_anchor=(0.4, -0.1), ncol=2, fontsize=8.5)
    title(fig, "What is hidden along with the test soil",
          f"Same {HOLDOUT_N:,} soils (those with depth and bulk density). Chance is 8 % "
          "exact, 25 % group.\nHiding the soil's whole data source costs about 9 points. "
          "Depth and bulk density add 4 points\nfor a new depth at a known site, "
          "under 2 (not significant) beyond it.")
    save(fig, "v6_holdout_levels")


def make_all(d, group, order):
    style()
    d = d.copy()
    rows = fig_limits(d, group)
    fig_confusion(d, group)
    fig_recall(d)
    fig_ks_scatter(d)
    fig_ks_cdf(d)
    fig_holdout()
    for r in rows:
        print(f"  {r['label']}: chance {r['chance']:.1f} %, nearest neighbour "
              f"{r['nn']:.1f} %, tool {r['base']:.1f} %, + depth & BD "
              f"{r['cov']:.1f} %, ceiling {r['hi']:.1f} %")
