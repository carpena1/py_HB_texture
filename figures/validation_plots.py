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
S3, S4, S5 = "#12968a", "#7b5bd6", "#c2255c"   # slots 3-5, for the methods

# The prediction methods compared in v1, v3 and v5 (cache column, label,
# colour). "heads" is the fixed-head predictor set (swcc_texture.HEADS_CM),
# the tool's default since 2026-09; the other arms use the vG parameters.
METHODS = [("pred_knn", "neighbour vote alone (vG)", S3),
           ("pred_base", "vG parameters, curve alone (former default)", S1),
           ("pred_cov", "vG parameters + depth & bulk density", S2),
           ("pred_heads_cov", "fixed heads + depth & BD (the default)", S4),
           ("pred_ens_cov", "vG ensemble of both opinions + depth & BD", S5)]
SEQ = LinearSegmentedColormap.from_list(
    "seq", [SURFACE, "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf",
            "#184f95", "#0d366b"])

# Classes ordered by texture group, so the groups form contiguous blocks.
GROUPED = ["sand", "loamy sand",
           "sandy loam", "loam", "sandy clay loam", "clay loam",
           "silt", "silt loam", "silty clay loam",
           "sandy clay", "silty clay", "clay"]

# Three methods at each hold-out level, each on its own script's targets, all
# with depth and bulk density: the vG parameters and the fixed heads from
# verify_alt_predictors.py (1,681 targets carrying both), and the ensemble
# from verify_blend.py (1,733 targets). (name, then exact/group per method.)
HOLDOUT_N = 1681
HOLDOUT = [("Soil alone\n(the standard)", 44.6, 66.8, 45.7, 68.8, 46.6, 67.9),
           ("Whole profile", 42.8, 66.9, 45.7, 68.5, 45.0, 67.6),
           ("Whole data source", 29.0, 57.0, 31.5, 58.8, 29.4, 56.9)]


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
    ax.set_xlabel("share of soils predicted correctly -- equal numbers per "
                  "class, so this is mean recall over the classes and takes "
                  "no account of precision", fontsize=8.5, color=INK2)
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
    for t, p in zip(d.texture_class, d.pred_heads_cov):
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
          "leave-one-soil-out, fixed heads (the default) with depth and bulk "
          "density.\nDiagonal = correct. "
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
    ax.set_xlabel("recall: share of the class that is found. A method can "
                  "raise this simply by predicting the class more often -- "
                  "v8 scores the same runs by F1", fontsize=8.5, color=INK2)
    h = [plt.Line2D([], [], marker="o", ls="", ms=9, color=c, mec=SURFACE,
                    mew=1.8) for _, _, c in METHODS]
    ax.legend(h, [lab for _, lab, _ in METHODS], loc="upper center",
              bbox_to_anchor=(0.45, -0.11), ncol=2, fontsize=8.5)
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
    ok = d.ksat_cmh.gt(0) & d.ks_heads_med.gt(0)
    sets = [("All soils with a measured Ks", d[ok])]
    fig, axes = plt.subplots(1, 1, figsize=(6.4, 5.6), squeeze=False)
    axes = axes[0]
    fig.subplots_adjust(top=0.8, bottom=0.12, left=0.14, right=0.84)
    lim = (-3.5, 2.5)
    xx = np.array(lim)
    for ax, (name, s) in zip(axes, sets):
        x, y = np.log10(s.ksat_cmh), np.log10(s.ks_heads_med)
        ax.fill_between(xx, xx - 1, xx + 1, color=GRID, alpha=0.55, lw=0, zorder=0)
        ax.fill_between(xx, xx - np.log10(2), xx + np.log10(2), color=GRID,
                        lw=0, zorder=0)
        ax.plot(xx, xx, color=INK2, lw=1, zorder=1)
        hb = ax.hexbin(x, y, gridsize=26, extent=lim + lim, mincnt=1, cmap=SEQ,
                       vmin=0, zorder=2, linewidths=0.2, edgecolors=SURFACE)
        f, w10, rho, n = ks_stats(s.ksat_cmh.to_numpy(), s.ks_heads_med.to_numpy())
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
          "Leave-one-soil-out, fixed heads (the default). Line = perfect "
          "prediction; darker band =\n"
          "within a factor of 2, lighter band = within a factor of 10.")
    save(fig, "v4_ks_predicted_vs_measured")


def fig_ks_cdf(d):
    """Share of soils whose Ks is within a factor F of the measurement."""
    cols = [("vG, curve only", "ks_med", S1),
            ("vG + bulk density", "ks_cov_med", S2),
            ("fixed heads (the default)", "ks_heads_med", S3),
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
    """Accuracy at each hold-out level for the three leading methods, all
    given depth and bulk density."""
    fig, ax = plt.subplots(figsize=(7.8, 5.6))
    fig.subplots_adjust(top=0.78, bottom=0.24, left=0.22, right=0.95)
    ys = np.arange(len(HOLDOUT))[::-1] * 1.55
    h = 0.17
    # exact class first (solid), then texture group (pale), in method order
    bars = ((S2, 1.0), (S4, 1.0), (S5, 1.0), (S2, 0.4), (S4, 0.4), (S5, 0.4))
    for y, (name, *v) in zip(ys, HOLDOUT):
        order = (v[0], v[2], v[4], v[1], v[3], v[5])
        for i, (val, (c, al)) in enumerate(zip(order, bars)):
            dy = (2.5 - i) * (h + 0.02)
            ax.barh(y + dy, val, height=h, color=c, alpha=al, zorder=2)
            ax.text(val + 1, y + dy, f"{val:.1f} %", va="center", fontsize=8,
                    color=INK)
    ax.set_yticks(ys, [n for n, *_ in HOLDOUT], fontsize=9.5)
    ax.set_xlim(0, 100)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0f} %")
    hgrid(ax, [0, 25, 50, 75, 100])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("share of soils predicted correctly -- equal numbers per "
                  "class, so this is mean recall and ignores precision "
                  "(v11 scores the same runs by F1)", fontsize=8.5,
                  color=INK2)
    hd = [Rectangle((0, 0), 1, 1, color=c, alpha=al) for c, al in bars]
    ax.legend(hd, ["exact, vG parameters", "exact, fixed heads",
                   "exact, ensemble", "group, vG parameters",
                   "group, fixed heads", "group, ensemble"],
              loc="upper center", bbox_to_anchor=(0.4, -0.12), ncol=3,
              fontsize=8.5)
    title(fig, "What is hidden along with the test soil, for the three leading methods",
          "About 1,700 soils, all given depth and bulk density. Chance is 8 % "
          "exact, 25 % group. Hiding the soil's whole data source\ncosts 14 to "
          "16 points. The fixed heads and the ensemble both beat the vG "
          "parameters when the soil's own source is in the\nreference; for a "
          "source the reference has never seen, only the fixed heads keep a "
          "clear margin.")
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


def per_class(d, col, c):
    """recall, precision, F1 (%) of one method on one class."""
    t = d.texture_class.to_numpy()
    p = d[col].to_numpy()
    tp = int(((p == c) & (t == c)).sum())
    n_true, n_pred = int((t == c).sum()), int((p == c).sum())
    rec = 100 * tp / n_true if n_true else np.nan
    prec = 100 * tp / n_pred if n_pred else np.nan
    f1 = 200 * tp / (n_true + n_pred) if n_true + n_pred else np.nan
    return rec, prec, f1


def macro_f1(d, col):
    return float(np.mean([per_class(d, col, c)[2] for c in GROUPED]))


def fig_f1(d):
    """Per-class F1 for each method, with the macro-F1 summary beside it.

    v3 shows recall alone, which a method can raise simply by predicting a
    class more often. F1 pairs it with precision, so a class only scores when
    the predictions of it are also right.
    """
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 6.2),
                             gridspec_kw={"width_ratios": [2.4, 1]})
    fig.subplots_adjust(top=0.82, bottom=0.2, left=0.17, right=0.98,
                        wspace=0.08)
    ax, ax2 = axes
    ys = np.arange(12)[::-1]
    for y, c in zip(ys, GROUPED):
        vals = [per_class(d, col, c)[2] for col, _, _ in METHODS]
        ax.plot([min(vals), max(vals)], [y, y], color=AXIS, lw=2, zorder=1)
        for v, (_, _, cc) in zip(vals, METHODS):
            ax.plot(v, y, "o", ms=8, color=cc, mec=SURFACE, mew=1.8, zorder=3)
    for b in (9.5, 5.5, 2.5):
        ax.axhline(b, color=GRID, lw=0.8)
    ax.set_yticks(ys, GROUPED, fontsize=9.5)
    ax.set_xlim(0, 90); ax.set_ylim(-0.7, 11.7)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0f}")
    hgrid(ax, [0, 25, 50, 75])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("F1 (harmonic mean of recall and precision)", fontsize=9,
                  color=INK2)

    names = [lab for _, lab, _ in METHODS]
    macro = [macro_f1(d, col) for col, _, _ in METHODS]
    order = np.argsort(macro)
    for k, i in enumerate(order):
        ax2.barh(k, macro[i], height=0.6, color=METHODS[i][2], zorder=2)
        ax2.text(macro[i] - 1, k, f"{macro[i]:.1f}", va="center", ha="right",
                 fontsize=9, color=SURFACE, fontweight="bold")
    ax2.set_yticks(range(len(METHODS)), [""] * len(METHODS))
    ax2.set_xlim(0, 55)
    ax2.set_title("macro-F1\n(the 12 classes, equally weighted)", fontsize=9,
                  color=INK2, pad=6)
    hgrid(ax2, [0, 25, 50])
    ax2.spines["left"].set_visible(False)
    ax2.tick_params(axis="y", length=0)

    h = [plt.Line2D([], [], marker="o", ls="", ms=9, color=c, mec=SURFACE,
                    mew=1.8) for _, _, c in METHODS]
    ax.legend(h, names, loc="upper center", bbox_to_anchor=(0.72, -0.1),
              ncol=2, fontsize=8.5)
    title(fig, "The same comparison as v3, scored by F1 instead of recall",
          "F1 asks whether the predictions of a class are also right. The "
          "neighbour vote's lead on the two rarest classes survives this\n"
          "(silt, sandy clay), but its lead elsewhere does not: it reaches "
          "rare classes by predicting them more often.")
    save(fig, "v8_f1_by_class")


def fig_precision_recall(d):
    """Why recall alone is not enough: the two move against each other."""
    show = [("pred_knn", "neighbour vote alone", S3),
            ("pred_cov", "vG parameters + depth & bulk density", S2),
            ("pred_ens_cov", "vG ensemble of both opinions", S5)]
    rare = ["silt", "sandy clay", "silty clay", "clay loam"]
    fig, ax = plt.subplots(figsize=(7.4, 6.2))
    fig.subplots_adjust(top=0.78, bottom=0.16, left=0.12, right=0.97)
    for f1 in (20, 40, 60, 80):                      # F1 contours
        r = np.linspace(f1 / 2 + 0.5, 100, 200)
        p = f1 * r / (2 * r - f1)
        ax.plot(r, p, color=GRID, lw=0.9, zorder=0)
        ax.text(99, min(f1 * 99 / (2 * 99 - f1), 99), f" F1 {f1}", fontsize=7.5,
                color=MUTED, va="center", ha="right", zorder=0)
    for c in GROUPED:
        pts = [per_class(d, col, c) for col, _, _ in show]
        xs = [q[0] for q in pts]
        ys = [q[1] for q in pts]
        ax.plot(xs, ys, color=AXIS, lw=0.8, zorder=1)
        for (rec, prec, _), (_, _, cc) in zip(pts, show):
            ax.plot(rec, prec, "o", ms=7.5, color=cc, mec=SURFACE, mew=1.5,
                    zorder=3)
        if c in rare:
            # label beside the highest-recall end, so the lines stay clear
            i = int(np.argmax(xs))
            ax.text(xs[i] + 1.8, ys[i] + 1.8, c, fontsize=8.5, color=INK,
                    ha="left", fontweight="bold", zorder=4)
    ax.plot([0, 100], [0, 100], color=GRID, lw=0.8, ls=":", zorder=0)
    ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    ax.set_xlabel("recall: share of the class that is found")
    ax.set_ylabel("precision: share of those predictions that are right")
    for t in (25, 50, 75):
        ax.axhline(t, color=GRID, lw=0.8, zorder=0)
    ax.set_xticks([0, 25, 50, 75, 100])
    h = [plt.Line2D([], [], marker="o", ls="", ms=9, color=c, mec=SURFACE,
                    mew=1.8) for _, _, c in show]
    ax.legend(h, [lab for _, lab, _ in show], loc="upper center",
              bbox_to_anchor=(0.5, -0.1), ncol=1, fontsize=9)
    title(fig, "Recall against precision: one line per texture class",
          "Each line joins the three methods for one class, leave-one-soil-out. "
          "Grey curves are equal F1. Moving right along a\nline buys recall "
          "with precision; only moving away from the origin is a real gain. "
          "The four labelled classes are the ones\nwhere the methods disagree "
          "most.")
    save(fig, "v9_precision_recall")


# verify_blend.py, both predictor sets and both members on the same targets
# in each design: (label, colour, {design: (exact, group, macro-F1)}).
BLEND = [
    ("vG parameters, classifier (former default)", S2,
     {"layer": (44.5, 66.6, 43.6), "profile": (43.6, 66.5, 42.1),
      "source": (29.0, 57.0, 27.0)}),
    ("vG parameters, neighbour vote", S3,
     {"layer": (42.5, 65.0, 42.2), "profile": (40.3, 63.1, 40.0),
      "source": (27.8, 54.5, 27.3)}),
    ("vG parameters, ensemble", S5,
     {"layer": (46.7, 67.8, 46.1), "profile": (45.2, 67.7, 44.7),
      "source": (29.9, 57.0, 28.7)}),
    ("fixed heads, classifier (the default)", S4,
     {"layer": (45.8, 68.6, 45.0), "profile": (45.8, 68.3, 45.0),
      "source": (31.4, 58.7, 29.7)}),
    ("fixed heads, ensemble", S1,
     {"layer": (47.2, 68.6, 46.6), "profile": (45.0, 68.1, 44.3),
      "source": (31.2, 58.7, 29.7)}),
]
DESIGNS = [("layer", "Soil alone\n(the standard)"),
           ("profile", "Whole profile"), ("source", "Whole data source")]


def macro_f1_group(d, col, group):
    """Macro-F1 over the four texture groups."""
    t = np.array([group[x] for x in d.texture_class])
    p = np.array([group[x] for x in d[col]])
    f = []
    for g in sorted(set(t)):
        tp = int(((p == g) & (t == g)).sum())
        den = int((p == g).sum() + (t == g).sum())
        f.append(200 * tp / den if den else 0.0)
    return float(np.mean(f))


def fig_macro_f1(d, group):
    """v1's layout, scored by macro-F1 instead of accuracy."""
    rows = []
    for label, f in (("Exact USDA class", lambda col: macro_f1(d, col)),
                     ("Texture group", lambda col: macro_f1_group(d, col,
                                                                  group))):
        rows.append(dict(label=label, nn=f("nn_class"),
                         **{col: f(col) for col, _, _ in METHODS}))

    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    fig.subplots_adjust(top=0.72, bottom=0.3, left=0.2, right=0.97)
    for y, r in zip((1, 0), rows):
        ax.plot([r["nn"]] * 2, [y - 0.34, y + 0.34], color=MUTED, lw=1.4)
        ax.text(r["nn"], y - 0.44, f"1 neighbour {r['nn']:.0f}", fontsize=8.5,
                color=INK2, ha="center", va="top")
        for i, (col, _, c) in enumerate(METHODS):
            dy = (i - 1.5) * 0.17
            ax.plot(r[col], y + dy, "o", ms=9, color=c, mec=SURFACE, mew=1.8,
                    zorder=4)
            ax.text(r[col] - 0.8, y + dy, f"{r[col]:.0f}", fontsize=8.5,
                    color=INK, ha="right", va="center")
    ax.set_yticks([1, 0], [r["label"] for r in rows], fontsize=10)
    ax.set_xlim(0, 100); ax.set_ylim(-0.75, 1.75)
    hgrid(ax, [0, 25, 50, 75, 100])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("macro-F1: recall and precision combined, classes weighted "
                  "equally. No ceiling is drawn -- the Cover & Hart bound "
                  "applies to the error rate, not to F1", fontsize=8.5,
                  color=INK2)
    h = [plt.Line2D([], [], marker="o", ls="", ms=9, color=c, mec=SURFACE,
                    mew=1.8) for _, _, c in METHODS]
    ax.legend(h, [lab for _, lab, _ in METHODS], loc="upper center",
              bbox_to_anchor=(0.42, -0.24), ncol=2, fontsize=8.5,
              handletextpad=0.4, columnspacing=1.2)
    title(fig, "The same methods as v1, scored by macro-F1",
          f"Leave-one-soil-out, {len(d):,} soils. v1 scores these runs by "
          "accuracy, which on balanced targets is mean recall; this one also\n"
          "charges for over-predicting a class. The ranking is the same, but "
          "the neighbour vote loses more ground than v1 suggests.")
    save(fig, "v10_macro_f1")


def fig_holdout_f1():
    """Macro-F1 at each hold-out level, both predictor sets, both members."""
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    fig.subplots_adjust(top=0.74, bottom=0.24, left=0.2, right=0.97)
    ys = np.arange(len(DESIGNS))[::-1] * 1.4
    h = 0.2
    for y, (key, _) in zip(ys, DESIGNS):
        for i, (lab, c, vals) in enumerate(BLEND):
            v = vals[key][2]
            dy = (2 - i) * (h + 0.02)
            ax.barh(y + dy, v, height=h, color=c, zorder=2)
            ax.text(v + 0.5, y + dy, f"{v:.1f}", va="center", fontsize=8.5,
                    color=INK)
    ax.set_yticks(ys, [lab for _, lab in DESIGNS], fontsize=9.5)
    ax.set_xlim(0, 55)
    hgrid(ax, [0, 10, 20, 30, 40, 50])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("macro-F1 (the 12 classes weighted equally)", fontsize=9,
                  color=INK2)
    hd = [Rectangle((0, 0), 1, 1, color=c) for _, c, _ in BLEND]
    ax.legend(hd, [lab for lab, _, _ in BLEND], loc="upper center",
              bbox_to_anchor=(0.42, -0.14), ncol=2, fontsize=8.5)
    title(fig, "What is hidden along with the test soil, scored by macro-F1",
          "1,733 soils per design, all given depth and bulk density "
          "(verify_blend.py). The ensemble lifts the vG parameters at every "
          "level,\nbut adds nothing once the fixed heads are used: the "
          "fixed-head classifier already captures what the neighbours "
          "contributed.")
    save(fig, "v11_f1_holdout_levels")


CM_METHODS = [("pred_knn", "neighbour vote alone"),
              ("pred_cov", "vG parameters + depth & BD"),
              ("pred_heads_cov", "fixed heads + depth & BD (the default)"),
              ("pred_ens_cov", "vG ensemble of both opinions")]


def confusion(d, col):
    """Row-normalised confusion matrix (% of each true class)."""
    idx = {c: i for i, c in enumerate(GROUPED)}
    m = np.zeros((12, 12))
    for t, p in zip(d.texture_class, d[col]):
        m[idx[t], idx[p]] += 1
    return m / np.maximum(m.sum(axis=1, keepdims=True), 1) * 100


def fig_confusion_methods(d):
    """One confusion matrix per method, on identical soils."""
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 11.4))
    fig.subplots_adjust(top=0.84, bottom=0.08, left=0.14, right=0.88,
                        hspace=0.5, wspace=0.18)
    for ax, (col, name) in zip(axes.ravel(), CM_METHODS):
        pct = confusion(d, col)
        im = ax.imshow(pct, cmap=SEQ, vmin=0, vmax=80)
        for i in range(12):
            v = pct[i, i]
            ax.text(i, i, f"{v:.0f}", ha="center", va="center", fontsize=7.5,
                    color="white" if v >= 45 else INK, fontweight="bold")
        for b in (1.5, 5.5, 8.5):
            ax.axhline(b, color=INK2, lw=0.7)
            ax.axvline(b, color=INK2, lw=0.7)
        ex = np.mean(d[col] == d.texture_class) * 100
        ax.set_title(f"{name}\nexact {ex:.1f} %, macro-F1 "
                     f"{macro_f1(d, col):.1f}", fontsize=9.5, color=INK,
                     pad=6)
        ax.set_xticks(range(12), GROUPED, rotation=90, fontsize=6.5)
        ax.set_yticks(range(12), GROUPED, fontsize=6.5)
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
    cax = fig.add_axes([0.9, 0.3, 0.018, 0.35])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("% of the true class", color=INK2)
    cb.outline.set_visible(False)
    cb.ax.tick_params(color=MUTED, labelcolor=INK2)
    fig.text(0.5, 0.055, "predicted class", ha="center", color=INK2,
             fontsize=9.5)
    fig.text(0.055, 0.46, "true class", rotation=90, va="center", color=INK2,
             fontsize=9.5)
    title(fig, "Where each method's errors go",
          f"Leave-one-soil-out, the same {len(d):,} soils in every panel, "
          "rows normalised to the true class. The diagonal is labelled; "
          "lines\nseparate the four texture groups. The off-diagonal weight "
          "sits next to the diagonal in every panel: the confusions are with "
          "neighbouring\ntextures, not with distant ones.")
    save(fig, "v12_confusion_by_method")


def fig_confusion_diff(d):
    """What changes when the fixed heads replace the vG parameters."""
    diff = confusion(d, "pred_heads_cov") - confusion(d, "pred_cov")
    lim = float(np.abs(diff).max())
    div = LinearSegmentedColormap.from_list(
        "div", ["#b8460f", "#e08a5a", "#f3d3bd", SURFACE, "#bcd6f3",
                "#5a8ed0", "#12447e"])
    fig, ax = plt.subplots(figsize=(8.2, 7.6))
    fig.subplots_adjust(top=0.84, bottom=0.2, left=0.2, right=0.93)
    im = ax.imshow(diff, cmap=div, vmin=-lim, vmax=lim)
    for i in range(12):
        for j in range(12):
            v = diff[i, j]
            if abs(v) >= 4:
                ax.text(j, i, f"{v:+.0f}", ha="center", va="center",
                        fontsize=7.5, color=INK,
                        fontweight="bold" if i == j else "normal")
    for b in (1.5, 5.5, 8.5):
        ax.axhline(b, color=INK2, lw=0.8)
        ax.axvline(b, color=INK2, lw=0.8)
    ax.set_xticks(range(12), GROUPED, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(12), GROUPED, fontsize=9)
    ax.set_xlabel("predicted class", color=INK2)
    ax.set_ylabel("true class", color=INK2)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("change in % of the true class", color=INK2)
    cb.outline.set_visible(False)
    cb.ax.tick_params(color=MUTED, labelcolor=INK2)
    title(fig, "What the fixed heads change, against the vG parameters",
          "The same two panels of v12, subtracted. Blue on the diagonal is a "
          "class predicted correctly more often; orange off it is a\n"
          "confusion the fixed heads make less. Cells under 4 points are "
          "unlabelled.")
    save(fig, "v13_confusion_change")


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
    fig_f1(d)
    fig_precision_recall(d)
    fig_macro_f1(d, group)
    fig_holdout_f1()
    fig_confusion_methods(d)
    fig_confusion_diff(d)
    for r in rows:
        methods = "  ".join(f"{lab}: {r[col]:.1f} %"
                            for col, lab, _ in METHODS)
        print(f"  {r['label']}: chance {r['chance']:.1f} %, nearest neighbour "
              f"{r['nn']:.1f} %, ceiling {r['hi']:.1f} %\n      {methods}")
