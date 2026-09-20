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
S3, S4 = "#12968a", "#7b5bd6"               # slots 3 and 4, for the methods

# The prediction methods compared in v1, v3 and v5 (cache column, label,
# colour). "heads" is the fixed-head predictor set of verify_alt_predictors.py,
# kept off by default in the tool.
METHODS = [("pred_knn", "neighbour vote alone", S3),
           ("pred_base", "hybrid: curve alone (the default)", S1),
           ("pred_cov", "hybrid + depth & bulk density", S2),
           ("pred_heads_cov", "fixed heads + depth & bulk density", S4)]
SEQ = LinearSegmentedColormap.from_list(
    "seq", [SURFACE, "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf",
            "#184f95", "#0d366b"])

# Classes ordered by texture group, so the groups form contiguous blocks.
GROUPED = ["sand", "loamy sand",
           "sandy loam", "loam", "sandy clay loam", "clay loam",
           "silt", "silt loam", "silty clay loam",
           "sandy clay", "silty clay", "clay"]

# verify_alt_predictors.py, all three designs on the same 1,681 targets that
# carry depth and bulk density, both methods through the shipped pipeline:
# (name, vG exact, vG group, heads exact, heads group), with the covariates.
HOLDOUT_N = 1681
HOLDOUT = [("Soil alone\n(the standard)", 44.6, 66.8, 45.7, 68.8),
           ("Whole profile", 42.8, 66.9, 45.7, 68.5),
           ("Whole data source", 29.0, 57.0, 31.5, 58.8)]


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
    """Every prediction method against chance and the Cover & Hart ceiling."""
    truth = d.texture_class.to_numpy()
    rows = []
    for label, n_cls, f in (("Exact USDA class", 12, lambda x: x),
                            ("Texture group", 4, lambda x: group[x])):
        t = np.array([f(x) for x in truth])
        acc = lambda col: np.mean(np.array([f(x) for x in d[col]]) == t) * 100
        nn = acc("nn_class")
        lo, hi = bayes_bracket(1 - nn / 100, n_cls)
        rows.append(dict(label=label, chance=100 / n_cls, nn=nn, hi=hi * 100,
                         **{col: acc(col) for col, _, _ in METHODS}))

    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    fig.subplots_adjust(top=0.72, bottom=0.26, left=0.2, right=0.97)
    for y, r in zip((1, 0), rows):
        ax.add_patch(Rectangle((r["nn"], y - 0.3), r["hi"] - r["nn"], 0.6,
                               color=GRID, zorder=1, lw=0))
        ax.plot([r["hi"]] * 2, [y - 0.34, y + 0.34], color=INK, lw=1.4, zorder=2)
        ax.text(r["hi"] + 0.8, y + 0.26, f"ceiling {r['hi']:.0f} %",
                fontsize=9, color=INK, va="center")
        ax.plot([r["chance"]] * 2, [y - 0.34, y + 0.34], color=MUTED, lw=1.4)
        ax.text(r["chance"], y - 0.44, f"chance {r['chance']:.0f} %",
                fontsize=8.5, color=INK2, ha="center", va="top")
        # one row of markers per method, stacked so labels never overprint
        for i, (col, _, c) in enumerate(METHODS):
            dy = (i - 1.5) * 0.17
            ax.plot(r[col], y + dy, "o", ms=9, color=c, mec=SURFACE, mew=1.8,
                    zorder=4)
            ax.text(r[col] - 1.2, y + dy, f"{r[col]:.0f}", fontsize=8.5,
                    color=INK, ha="right", va="center")
    ax.set_yticks([1, 0], [r["label"] for r in rows], fontsize=10)
    ax.set_xlim(0, 100); ax.set_ylim(-0.75, 1.75)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0f} %")
    hgrid(ax, [0, 25, 50, 75, 100])
    ax.spines["left"].set_visible(False)
    h = [plt.Line2D([], [], marker="o", ls="", ms=9, color=c, mec=SURFACE,
                    mew=1.8) for _, _, c in METHODS]
    h.append(Rectangle((0, 0), 1, 1, color=GRID))
    ax.legend(h, [lab for _, lab, _ in METHODS] +
              ["achievable range (nearest neighbour to ceiling)"],
              loc="upper center", bbox_to_anchor=(0.42, -0.2), ncol=2,
              fontsize=8.5, handletextpad=0.4, columnspacing=1.2)
    title(fig, "Every method tried, against the best achievable accuracy",
          f"Leave-one-soil-out, {len(d):,} soils (equal numbers per class). "
          "The grey band runs from plain nearest-neighbour matching to the\n"
          "Cover & Hart (1967) bound, which is the same for all four methods: "
          "they re-encode one fitted curve, and cannot add information.")
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
    """Per-class recall for each prediction method."""
    fig, ax = plt.subplots(figsize=(7.8, 6.2))
    fig.subplots_adjust(top=0.85, bottom=0.18, left=0.2, right=0.9)
    ys = np.arange(12)[::-1]
    for y, c in zip(ys, GROUPED):
        m = d.texture_class == c
        vals = [np.mean(d[col][m] == c) * 100 for col, _, _ in METHODS]
        ax.plot([min(vals), max(vals)], [y, y], color=AXIS, lw=2, zorder=1)
        for v, (_, _, col_c) in zip(vals, METHODS):
            ax.plot(v, y, "o", ms=8, color=col_c, mec=SURFACE, mew=1.8, zorder=3)
        ax.text(101, y, f"n={int(m.sum())}", fontsize=8, color=MUTED,
                va="center")
    ax.axvline(100 / 12, color=MUTED, lw=1.2)
    ax.text(100 / 12 + 1, 11.55, "chance", fontsize=8.5, color=INK2)
    for b in (9.5, 5.5, 2.5):
        ax.axhline(b, color=GRID, lw=0.8)
    ax.set_yticks(ys, GROUPED, fontsize=9.5)
    ax.set_xlim(0, 100); ax.set_ylim(-0.7, 11.9)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0f} %")
    hgrid(ax, [0, 25, 50, 75, 100])
    ax.spines["left"].set_visible(False)
    h = [plt.Line2D([], [], marker="o", ls="", ms=9, color=c, mec=SURFACE,
                    mew=1.8) for _, _, c in METHODS]
    ax.legend(h, [lab for _, lab, _ in METHODS], loc="upper center",
              bbox_to_anchor=(0.45, -0.09), ncol=2, fontsize=8.5)
    title(fig, "Correct class, by true class, for each method",
          "Share of each class predicted exactly, leave-one-soil-out. "
          "The methods differ most where the class is rare:\nsilt, sandy clay "
          "and silty clay. Horizontal lines separate the texture groups.")
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
    cols = [("retention curve only", "ks_med", S1),
            ("+ bulk density", "ks_cov_med", S2),
            ("fixed heads", "ks_heads_med", S3),
            ("fixed heads + bulk density", "ks_heads_cov_med", S4)]
    ok = d.ksat_cmh.gt(0)
    for _, c, _ in cols:
        ok &= d[c].gt(0)
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    fig.subplots_adjust(top=0.78, bottom=0.14, left=0.1, right=0.8)
    F = np.logspace(0, 3, 300)
    rows = []
    for name, col, c in cols:
        s_ = d[ok]
        e = np.abs(np.log10(s_[col]) - np.log10(s_.ksat_cmh)).to_numpy()
        ax.plot(F, [np.mean(e <= np.log10(f)) * 100 for f in F], color=c, lw=2,
                solid_capstyle="round", label=name)
        vals = [np.mean(e <= np.log10(f)) * 100 for f in (2, 10)]
        for f, v in zip((2, 10), vals):
            ax.plot(f, v, "o", ms=7, color=c, mec=SURFACE, mew=2, zorder=3)
        rows.append(f"{name:<28s}{vals[0]:4.0f} %{vals[1]:7.0f} %")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.26), fontsize=8.5)
    ax.text(0.97, 0.03, f"{'':<28s}{'2x':>6s}{'10x':>9s}\n" + "\n".join(rows),
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            family="monospace", color=INK, linespacing=1.5)
    ax.set_xscale("log")
    ax.set_xlim(1, 1000); ax.set_ylim(0, 100)
    ax.set_xticks([1, 2, 5, 10, 100, 1000], ["1", "2", "5", "10", "100", "1000"])
    ax.yaxis.set_major_formatter(lambda v, _: f"{v:.0f} %")
    for t in (25, 50, 75, 100):
        ax.axhline(t, color=GRID, lw=0.8, zorder=0)
    ax.set_xlabel("within a factor of ... of the measured Ks")
    ax.set_ylabel("share of soils")
    title(fig, "How close the Ks estimate is, by method",
          "Share of soils whose predicted Ks lies within a given factor of the "
          "measurement, leave-one-soil-out.\nBulk density helps; the fixed-head "
          "predictors do not change Ks when the soil's own source is in the "
          "reference.")
    save(fig, "v5_ks_within_factor")


def fig_holdout():
    """Accuracy at each hold-out level, vG parameters against fixed heads,
    both with depth and bulk density supplied (verify_alt_predictors.py)."""
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    fig.subplots_adjust(top=0.76, bottom=0.22, left=0.22, right=0.95)
    ys = np.arange(len(HOLDOUT))[::-1] * 1.25
    h = 0.24
    bars = ((S2, 1.0), (S4, 1.0), (S2, 0.45), (S4, 0.45))
    for y, (name, *vals) in zip(ys, HOLDOUT):
        # vals = vG exact, vG group, heads exact, heads group
        for i, (v, (c, al)) in enumerate(zip((vals[0], vals[2], vals[1],
                                              vals[3]), bars)):
            dy = (1.5 - i) * (h + 0.02)
            ax.barh(y + dy, v, height=h, color=c, alpha=al, zorder=2)
            ax.text(v + 1, y + dy, f"{v:.1f} %", va="center", fontsize=8.5,
                    color=INK)
    ax.set_yticks(ys, [n for n, *_ in HOLDOUT], fontsize=9.5)
    ax.set_xlim(0, 100)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0f} %")
    hgrid(ax, [0, 25, 50, 75, 100])
    ax.spines["left"].set_visible(False)
    hd = [Rectangle((0, 0), 1, 1, color=c, alpha=al) for c, al in bars]
    ax.legend(hd, ["exact class, vG parameters", "exact class, fixed heads",
                   "texture group, vG parameters", "texture group, fixed heads"],
              loc="upper center", bbox_to_anchor=(0.4, -0.12), ncol=2,
              fontsize=8.5)
    title(fig, "What is hidden along with the test soil, for both predictor sets",
          f"Same {HOLDOUT_N:,} soils (those with depth and bulk density), both "
          "given depth and bulk density. Chance is 8 % exact, 25 % group.\n"
          "Hiding the soil's whole data source costs 14 to 16 points. The fixed "
          "heads lead at every level, by 1.1 to 2.9 points, most where it "
          "matters\nmost: a data source the reference has never seen.")
    save(fig, "v6_holdout_levels")


def fig_fractions(d):
    """Particle fractions: the neighbour mean against the regressors, in the
    reference and on outside sources."""
    parts = ("sand", "silt", "clay")
    inref = [[np.nanmean(np.abs(d[f"{pre}_{c}"] - d[c])) for c in parts]
             for pre in ("fknn", "fgbm_cov")]
    # verify_fractions.py --cov: five external sets, each held out of the
    # reference and predicted from its measured points.
    outside = [[15.27, 11.03, 10.90], [15.94, 10.29, 14.31]]
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.2), sharey=True)
    fig.subplots_adjust(top=0.72, bottom=0.22, left=0.09, right=0.98,
                        wspace=0.12)
    for ax, vals, name, n in (
            (axes[0], inref, "soils from the reference itself", len(d)),
            (axes[1], outside, "five outside sources, each held out", 378)):
        x = np.arange(3)
        for k, (v, c, lab) in enumerate(zip(vals, (S1, S2),
                                            ("neighbour mean (the default)",
                                             "sand & clay regressors"))):
            ax.bar(x + (k - 0.5) * 0.36, v, width=0.34, color=c, zorder=2,
                   label=lab)
            for xi, vi in zip(x + (k - 0.5) * 0.36, v):
                ax.text(xi, vi + 0.3, f"{vi:.1f}", ha="center", fontsize=8.5,
                        color=INK)
        ax.set_xticks(x, parts, fontsize=10)
        ax.set_title(f"{name}  (n = {n:,})", fontsize=9.5, color=INK2, pad=6)
        for t in (5, 10, 15):
            ax.axhline(t, color=GRID, lw=0.8, zorder=0)
        ax.set_ylim(0, 18)
        ax.tick_params(axis="x", length=0)
    axes[0].set_ylabel("mean absolute error (points)", color=INK2)
    axes[0].legend(loc="upper center", bbox_to_anchor=(1.0, -0.12), ncol=2,
                   fontsize=9)
    title(fig, "Particle fractions: why the regressors were not adopted",
          "Predicting sand and clay directly, with silt by difference, beats "
          "the neighbour mean on soils drawn from the reference\n(left) and "
          "loses on a source the reference has never seen (right), where "
          "boosting shrinks extreme clays towards the mean.")
    save(fig, "v7_fraction_methods")


def make_all(d, group, order):
    style()
    d = d.copy()
    rows = fig_limits(d, group)
    fig_confusion(d, group)
    fig_recall(d)
    fig_ks_scatter(d)
    fig_ks_cdf(d)
    fig_holdout()
    fig_fractions(d)
    for r in rows:
        methods = "  ".join(f"{lab}: {r[col]:.1f} %"
                            for col, lab, _ in METHODS)
        print(f"  {r['label']}: chance {r['chance']:.1f} %, nearest neighbour "
              f"{r['nn']:.1f} %, ceiling {r['hi']:.1f} %\n      {methods}")
