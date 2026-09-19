"""Per-soil report for an external test set: predicted against reported texture,
fractions and Ks, as tables and figures.

Every soil is predicted from its measured retention points with the default
reference (the set itself left out), as the command line would, under four
covariate options: the curve only, + bulk density, + organic carbon, and
both. Bulk density reaches the classifier and the neighbour search; organic
carbon reaches the neighbour search only (fractions and Ks), the classifier
does not take it. For each output -- class, fractions, Ks -- the option that
scores best over the set is reported, and named in the table's covariates
columns. That choice is made on the same soils it is scored on, so the
reported numbers are an upper bound; all four options are in the stats table.

Outputs go to figures/external/<name>/ (git-ignored while the set is not
distributed): soils.csv, groups.csv, stats.csv, report.md and fig1-fig5.

Usage:  python figures/make_external_report.py babaeian_az [n_mc] [--replot]

--replot redraws from the saved predictions_all_options.csv.
"""

import ast
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "figures"))
os.chdir(ROOT)

import matplotlib                                              # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                # noqa: E402
from scipy.stats import spearmanr                              # noqa: E402

import ksat_metrics as km                                      # noqa: E402
import swcc_texture as st                                      # noqa: E402
import validation_plots as vp                                  # noqa: E402
from verify_common import GROUP, GROUP_ORDER                   # noqa: E402
from verify_external import DATASETS                           # noqa: E402

ARMS = {  # name: (bulk density, organic carbon)
    "none": (False, False),
    "bulk density": (True, False),
    "organic carbon": (False, True),
    "bulk density + organic carbon": (True, True),
}
# Sets whose reference Ks is empty are scored against another measured Ks,
# named here, for information only.
KS_TRUTH = {"willard": ("ksat_mpd_cmh",
                        "field permeameter Ks, dry season (not in the reference)")}
LABELS = {"babaeian_az": "Arizona soils", "willard": "Laikipia soils",
          "boorowa": "Boorowa Farm soils (CSIRO, NSW)"}
CURVE_NOTE = {"willard": "saturation at 95 % of porosity and the pressure-plate "
                         "points, pF 2.15-4.2",
              "boorowa": "10 cm to 15 bar: suction tables and pressure plates"}
MATCH_COLOR = {"exact": "#2a78d6", "same group": "#e8a33d",
               "wrong group": "#d6452a"}


def predict(name, n_mc):
    tg = st.load_reference_df(name).reset_index(drop=True)
    pts = DATASETS[name].measured_points()
    ref_df = st.load_reference_df()
    ref_df = ref_df[~ref_df.layer_id.isin(tg.layer_id)].reset_index(drop=True)
    clfs = {False: st.TextureGBM(df=ref_df),
            True: st.TextureGBM(df=ref_df, covariates=["bd"])}
    rows = []
    for arm, (use_bd, use_oc) in ARMS.items():
        ref = st.GshpReference(df=ref_df, use_bd=use_bd, use_om=use_oc)
        for _, r in tg.iterrows():
            h, th = pts[r.layer_id]
            e = st.estimate(h, th, ref=ref, clf=clfs[use_bd], n_mc=n_mc,
                            sample_type=r.sample_type,
                            bulk_density=r.bd if use_bd else None,
                            om=r.oc if use_oc else None)   # the reference matches on oc
            f, k = e["fractions"], e["ksat"]
            probs = e["class_probabilities"]
            rows.append(dict(
                arm=arm, layer_id=r.layer_id, pred_class=e["texture_class"],
                pred_prob=probs[e["texture_class"]], top2=list(probs)[:2],
                **{f"pred_{c}": f[c] for c in ("sand", "silt", "clay")},
                **{f"pred_{c}_p5": f["p5"][c] for c in ("sand", "silt", "clay")},
                **{f"pred_{c}_p95": f["p95"][c] for c in ("sand", "silt", "clay")},
                ks_med=k["median_cmh"], ks_p5=k["p5_cmh"], ks_p95=k["p95_cmh"]))
    return tg, pd.DataFrame(rows)


def arm_stats(tg, p):
    """One row of set statistics per covariate option."""
    out = []
    for arm, g in p.groupby("arm", sort=False):
        g = g.set_index("layer_id").loc[tg.layer_id]
        truth = tg.texture_class.to_numpy()
        pred = g.pred_class.to_numpy()
        row = dict(covariates=arm, n=len(tg),
                   exact_pct=np.mean(pred == truth) * 100,
                   group_pct=np.mean([GROUP[a] == GROUP[b]
                                      for a, b in zip(pred, truth)]) * 100,
                   top2_pct=np.mean([t in s for t, s in zip(truth, g.top2)]) * 100)
        # Soils an option cannot run (organic carbon missing) are left out of
        # its fraction statistics; fractions_n says how many were scored.
        has = np.isfinite(g.pred_sand.to_numpy(float))
        row["fractions_n"] = int(has.sum())
        for c in ("sand", "silt", "clay"):
            obs = tg[c].to_numpy()[has]
            err = g[f"pred_{c}"].to_numpy()[has] - obs
            ss = np.sum((obs - obs.mean()) ** 2)
            row.update({f"{c}_mae": np.mean(np.abs(err)),
                        f"{c}_bias": np.mean(err),
                        f"{c}_rmse": np.sqrt(np.mean(err ** 2)),
                        f"{c}_r2": 1 - np.sum(err ** 2) / ss,
                        f"{c}_cover_pct": np.mean(
                            (obs >= g[f"pred_{c}_p5"].to_numpy()[has])
                            & (obs <= g[f"pred_{c}_p95"].to_numpy()[has])) * 100})
        row["fractions_mae"] = np.mean([row[f"{c}_mae"] for c in ("sand", "silt", "clay")])
        s = km.ksat_scores(tg.ksat_cmh.to_numpy(float), g.ks_med.to_numpy(float),
                           g.ks_p5.to_numpy(float), g.ks_p95.to_numpy(float))
        e = np.abs(np.log10(g.ks_med.to_numpy(float)) - np.log10(tg.ksat_cmh.to_numpy(float)))
        row.update(ks_n=s["n"], ks_median_factor=10 ** np.nanmedian(e), ks_bias_log10=s["bias_log10"],
                   ks_rmse_log10=s["rmse_log10"], ks_within_2x_pct=s["within_2x"] * 100,
                   ks_within_5x_pct=s["within_5x"] * 100,
                   ks_within_10x_pct=s["within_10x"] * 100,
                   ks_cover_pct=s["coverage_p5_p95"] * 100, ks_spearman=s["spearman"])
        out.append(row)
    return pd.DataFrame(out)


def best_arms(stats):
    """Class: exact, then group, then top-2 accuracy. Fractions: mean absolute
    error, among the options that score every soil. Ks: median |log10 error|,
    then RMSE. Ties keep the fewer covariates."""
    s = stats.reset_index(drop=True)
    cls = s.sort_values(["exact_pct", "group_pct", "top2_pct"], ascending=False,
                        kind="stable").covariates.iloc[0]
    full = s[s.fractions_n == s.n]
    frac = full.sort_values("fractions_mae", kind="stable").covariates.iloc[0]
    ks = s.sort_values(["ks_median_factor", "ks_rmse_log10"],
                       kind="stable").covariates.iloc[0]
    return cls, frac, ks


def soil_table(tg, p, cls, frac, ks):
    pc = p[p.arm == cls].set_index("layer_id")
    pf = p[p.arm == frac].set_index("layer_id")
    pk = p[p.arm == ks].set_index("layer_id")
    rows = []
    for _, r in tg.iterrows():
        c, f, k = pc.loc[r.layer_id], pf.loc[r.layer_id], pk.loc[r.layer_id]
        match = ("exact" if c.pred_class == r.texture_class else
                 "same group" if GROUP[c.pred_class] == GROUP[r.texture_class]
                 else "wrong group")
        rows.append({
            "sample": r.layer_id.split("_", 1)[1],
            "reported class": r.texture_class, "predicted class": c.pred_class,
            "class probability": round(c.pred_prob, 2),
            "reported group": GROUP[r.texture_class],
            "predicted group": GROUP[c.pred_class], "match": match,
            "covariates (class)": cls,
            **{f"{x} reported": round(r[x], 1) for x in ("sand", "silt", "clay")},
            **{f"{x} predicted [p5-p95]":
               f"{f[f'pred_{x}']:.1f} [{f[f'pred_{x}_p5']:.0f}-{f[f'pred_{x}_p95']:.0f}]"
               for x in ("sand", "silt", "clay")},
            "covariates (fractions)": frac,
            "Ks measured (cm/h)": round(r.ksat_cmh, 3),
            "Ks predicted (cm/h)": round(k.ks_med, 3),
            "Ks p5-p95 (cm/h)": f"{k.ks_p5:.3g}-{k.ks_p95:.3g}",
            "Ks predicted/measured": round(k.ks_med / r.ksat_cmh, 2),
            "covariates (Ks)": ks,
            "bulk density (g/cm3)": r.bd, "organic matter (%)": round(r.om, 2),
            "sample type": r.sample_type,
            "_clay_order": r.clay})
    return (pd.DataFrame(rows).sort_values(["_clay_order"])
            .drop(columns="_clay_order").reset_index(drop=True))


def group_table(soils):
    conf = pd.crosstab(pd.Categorical(soils["reported group"], GROUP_ORDER),
                       pd.Categorical(soils["predicted group"], GROUP_ORDER),
                       dropna=False)
    conf.index.name, conf.columns.name = "reported group", "predicted group"
    rows = []
    for gname in GROUP_ORDER:
        n_rep = int(conf.loc[gname].sum())
        n_pred = int(conf[gname].sum())
        hit = int(conf.loc[gname, gname])
        rows.append({"group": gname, "reported": n_rep, "predicted": n_pred,
                     "correct": hit,
                     "recall %": round(hit / n_rep * 100, 1) if n_rep else np.nan,
                     "precision %": round(hit / n_pred * 100, 1) if n_pred else np.nan})
    summary = pd.DataFrame(rows)
    return conf, summary


# --------------------------------------------------------------------------
# USDA texture triangle (clay, silt, sand in %), polygons after the NRCS chart.
USDA_POLY = {
    "clay": [(100, 0, 0), (55, 0, 45), (40, 15, 45), (40, 40, 20), (60, 40, 0)],
    "silty clay": [(60, 40, 0), (40, 40, 20), (40, 60, 0)],
    "silty clay loam": [(40, 60, 0), (40, 40, 20), (27, 53, 20), (27, 73, 0)],
    "sandy clay": [(55, 0, 45), (35, 0, 65), (35, 20, 45), (40, 15, 45)],
    "clay loam": [(40, 40, 20), (40, 15, 45), (27, 28, 45), (27, 53, 20)],
    "sandy clay loam": [(35, 0, 65), (20, 0, 80), (20, 28, 52), (27, 28, 45),
                        (35, 20, 45)],
    "loam": [(27, 28, 45), (20, 28, 52), (7, 41, 52), (7, 50, 43), (27, 50, 23)],
    "silt loam": [(27, 50, 23), (7, 50, 43), (0, 50, 50), (0, 80, 20),
                  (12, 80, 8), (12, 88, 0), (27, 73, 0)],
    "silt": [(12, 88, 0), (12, 80, 8), (0, 80, 20), (0, 100, 0)],
    "sandy loam": [(15, 0, 85), (20, 0, 80), (20, 28, 52), (7, 41, 52),
                   (7, 50, 43), (0, 50, 50), (0, 30, 70)],
    "loamy sand": [(10, 0, 90), (15, 0, 85), (0, 30, 70), (0, 15, 85)],
    "sand": [(0, 0, 100), (10, 0, 90), (0, 15, 85)],
}


# Label positions (clay, silt, sand) for classes too small for their centroid.
LABEL_AT = {"sand": (2, 3, 95), "loamy sand": (5, 10, 85),
            "sandy loam": (12, 18, 70), "sandy clay loam": (27, 13, 60),
            "silt": (5, 88, 7), "clay": (70, 15, 15)}


def tern(clay, silt, sand):
    tot = np.asarray(clay) + np.asarray(silt) + np.asarray(sand)
    clay, silt = np.asarray(clay) / tot, np.asarray(silt) / tot
    return silt + clay / 2, clay * np.sqrt(3) / 2


def save(fig, out, name):
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=200,
                    bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print(f"  wrote {out}/{name}.png/.svg")


def fig_triangle(tg, pf, soils, out, frac_arm, cls_arm, label):
    fig, ax = plt.subplots(figsize=(7.6, 6.6))
    fig.subplots_adjust(top=0.88, bottom=0.02, left=0.02, right=0.98)
    for cname, poly in USDA_POLY.items():
        x, y = tern(*np.array(poly + poly[:1]).T)
        ax.plot(x, y, color=vp.AXIS, lw=0.8, zorder=1)
        cx, cy = tern(*LABEL_AT.get(cname, np.array(poly).mean(axis=0)))
        ax.text(cx, cy, cname, ha="center", va="center", fontsize=7.5,
                color=vp.MUTED, zorder=1)
    match = soils.set_index("sample")["match"]
    for _, r in tg.iterrows():
        s = r.layer_id.split("_", 1)[1]
        f = pf.loc[r.layer_id]
        x0, y0 = tern(r.clay, r.silt, r.sand)
        x1, y1 = tern(f.pred_clay, f.pred_silt, f.pred_sand)
        col = MATCH_COLOR[match[s]]
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0), zorder=2,
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=0.9,
                                    shrinkA=3, shrinkB=3, mutation_scale=8))
        ax.plot(x0, y0, "o", ms=6, color=col, zorder=3)
        ax.plot(x1, y1, "o", ms=6, mfc=vp.SURFACE, mec=col, mew=1.2, zorder=3)
        if len(tg) <= 30:    # too many to label legibly
            ax.text(x0 + 0.012, y0 + 0.008, s, fontsize=6.5, color=vp.INK2, zorder=4)
    for t in (20, 40, 60, 80):
        # tick labels: clay on the left edge, silt on the right, sand below
        x, y = tern(t, 0, 100 - t); ax.text(x - 0.02, y, f"{t}", ha="right", va="center", fontsize=7, color=vp.MUTED)
        x, y = tern(100 - t, t, 0); ax.text(x + 0.02, y, f"{t}", ha="left", va="center", fontsize=7, color=vp.MUTED)
        x, y = tern(0, 100 - t, t); ax.text(x, y - 0.03, f"{t}", ha="center", va="top", fontsize=7, color=vp.MUTED)
    ax.text(0.16, 0.47, "clay %", rotation=60, ha="center", color=vp.INK2, fontsize=9)
    ax.text(0.84, 0.47, "silt %", rotation=-60, ha="center", color=vp.INK2, fontsize=9)
    ax.text(0.5, -0.08, "sand %", ha="center", color=vp.INK2, fontsize=9)
    handles = [plt.Line2D([], [], marker="o", ls="", color=vp.INK2, label="reported"),
               plt.Line2D([], [], marker="o", ls="", mfc=vp.SURFACE, mec=vp.INK2,
                          label="predicted fractions")]
    handles += [plt.Line2D([], [], color=c, lw=2, label={
        "exact": "class predicted exactly", "same group": "right texture group only",
        "wrong group": "wrong group"}[m]) for m, c in MATCH_COLOR.items()]
    ax.legend(handles=handles, loc="upper right", fontsize=8.5, bbox_to_anchor=(1.0, 0.98))
    ax.set_xlim(-0.08, 1.08); ax.set_ylim(-0.11, 0.89)
    ax.set_aspect("equal"); ax.axis("off")
    vp.title(fig, f"{label}: reported and predicted texture",
             f"Arrow from the measured fractions to the predicted ones "
             f"(covariates: {frac_arm}).\nColour: whether the predicted class "
             f"(covariates: {cls_arm}) is exact, in the same group, or wrong.")
    save(fig, out, "fig1_texture_triangle")


def fig_confusion(soils, conf, out, cls_arm, label):
    order = [c for c in vp.GROUPED if c in set(soils["reported class"]) | set(soils["predicted class"])]
    m = pd.crosstab(pd.Categorical(soils["reported class"], order),
                    pd.Categorical(soils["predicted class"], order), dropna=False)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 5.4),
                                 gridspec_kw=dict(width_ratios=[1.45, 1]))
    fig.subplots_adjust(top=0.78, bottom=0.2, left=0.12, right=0.98, wspace=0.3)
    for ax, mat, labels, groups in ((a1, m.to_numpy(), order, True),
                                    (a2, conf.to_numpy(), GROUP_ORDER, False)):
        vmax = max(mat.max(), 1)
        ax.imshow(mat, cmap=vp.SEQ, vmin=0, vmax=vmax * 1.15)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                if mat[i, j]:
                    ax.text(j, i, int(mat[i, j]), ha="center", va="center", fontsize=9,
                            color="white" if mat[i, j] > vmax * 0.55 else vp.INK)
                if i == j:
                    ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False,
                                               ec=vp.INK2, lw=0.8))
        ax.set_xticks(range(len(labels)), labels, rotation=40, ha="right")
        ax.set_yticks(range(len(labels)), labels)
        ax.set_xlabel("predicted"); ax.set_ylabel("reported")
        ax.spines[:].set_visible(False)
        ax.tick_params(length=0)
    exact = np.trace(m.to_numpy()) / m.to_numpy().sum() * 100
    grp = np.trace(conf.to_numpy()) / conf.to_numpy().sum() * 100
    a1.set_title(f"USDA class — {exact:.0f} % exact", fontsize=10, color=vp.INK, loc="left")
    a2.set_title(f"Texture group — {grp:.0f} % right", fontsize=10, color=vp.INK, loc="left")
    vp.title(fig, f"{label}: reported against predicted class",
             f"Counts of soils; the diagonal is a correct prediction. "
             f"Covariates: {cls_arm}.")
    save(fig, out, "fig2_class_confusion")


def fig_fractions(tg, pf, out, frac_arm, label):
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.6))
    fig.subplots_adjust(top=0.78, bottom=0.14, left=0.06, right=0.98, wspace=0.28)
    for ax, c in zip(axes, ("sand", "silt", "clay")):
        f = pf.loc[tg.layer_id]
        x, y = tg[c].to_numpy(), f[f"pred_{c}"].to_numpy()
        lo, hi = f[f"pred_{c}_p5"].to_numpy(), f[f"pred_{c}_p95"].to_numpy()
        ax.plot([0, 100], [0, 100], color=vp.INK2, lw=1, zorder=1)
        ax.errorbar(x, y, yerr=[y - lo, hi - y], fmt="o", ms=5, color=vp.S1,
                    ecolor="#9ec5f4", elinewidth=1.2, zorder=2)
        err = y - x
        r2 = 1 - np.sum(err ** 2) / np.sum((x - x.mean()) ** 2)
        cover = np.mean((x >= lo) & (x <= hi)) * 100
        ax.text(0.03, 0.97, f"mean |error| {np.mean(np.abs(err)):.1f} pts\n"
                f"bias {np.mean(err):+.1f}   RMSE {np.sqrt(np.mean(err**2)):.1f}\n"
                f"R² {r2:.2f}   in p5–p95: {cover:.0f} %",
                transform=ax.transAxes, va="top", fontsize=8.5, color=vp.INK,
                linespacing=1.4, bbox=dict(facecolor=vp.SURFACE, edgecolor="none", pad=3))
        top = max(x.max(), hi.max()) * 1.08 + 2
        ax.set_xlim(0, top); ax.set_ylim(0, top); ax.set_aspect("equal")
        ax.set_xlabel(f"reported {c} (%)"); ax.set_ylabel(f"predicted {c} (%)")
    vp.title(fig, f"{label}: particle fractions",
             f"Predicted mean with its 5–95 % range against the reported value. "
             f"Covariates: {frac_arm}.")
    save(fig, out, "fig3_fractions")


def fig_ks(tg, pk, out, ks_arm, label, ks_label="measured Ks"):
    fig, ax = plt.subplots(figsize=(6.8, 6.2))
    fig.subplots_adjust(top=0.84, bottom=0.11, left=0.14, right=0.97)
    k = pk.loc[tg.layer_id]
    ok = (np.isfinite(tg.ksat_cmh.to_numpy(float))
          & np.isfinite(k.ks_med.to_numpy(float)))
    tg, k = tg[ok], k[ok]
    x, y = tg.ksat_cmh.to_numpy(float), k.ks_med.to_numpy(float)
    lo, hi = k.ks_p5.to_numpy(float), k.ks_p95.to_numpy(float)
    lim = (10 ** np.floor(np.log10(min(x.min(), lo.min()))),
           10 ** np.ceil(np.log10(max(x.max(), hi.max()))))
    xx = np.array(lim)
    ax.fill_between(xx, xx / 10, xx * 10, color=vp.GRID, alpha=0.55, lw=0, zorder=0)
    ax.fill_between(xx, xx / 2, xx * 2, color=vp.GRID, lw=0, zorder=0)
    ax.plot(xx, xx, color=vp.INK2, lw=1, zorder=1)
    gcol = dict(zip(GROUP_ORDER, ["#e8a33d", "#2a78d6", "#3aa37a", "#9b59b6"]))
    for gname in GROUP_ORDER:
        s = tg.texture_class.map(GROUP).to_numpy() == gname
        if s.any():
            ax.errorbar(x[s], y[s], yerr=[y[s] - lo[s], hi[s] - y[s]], fmt="o", ms=6,
                        color=gcol[gname], ecolor=gcol[gname], alpha=0.9,
                        elinewidth=0.8, capsize=0, label=gname, zorder=3)
    for xi, yi, lid in zip(x, y, tg.layer_id):
        if len(x) > 30:      # too many to label legibly
            break
        ax.text(xi * 1.1, yi * 1.05, lid.split("_", 1)[1], fontsize=6.5, color=vp.INK2, zorder=4)
    e = np.log10(y) - np.log10(x)
    ax.text(0.97, 0.03,
            f"n = {len(x)}\ntypical error ×{10 ** np.median(np.abs(e)):.1f}\n"
            f"bias {np.mean(e):+.2f} log10   RMSE {np.sqrt(np.mean(e**2)):.2f}\n"
            f"within 2× {np.mean(np.abs(e) <= np.log10(2))*100:.0f} %   "
            f"within 10× {np.mean(np.abs(e) <= 1)*100:.0f} %\n"
            f"in p5–p95: {np.mean((x >= lo) & (x <= hi))*100:.0f} %\n"
            f"rank correlation ρ = {spearmanr(x, y).statistic:.2f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5,
            color=vp.INK, linespacing=1.4,
            bbox=dict(facecolor=vp.SURFACE, edgecolor="none", alpha=0.9, pad=4))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
    ax.set_xlabel(f"{ks_label} (cm/h)"); ax.set_ylabel("predicted Ks (cm/h)")
    ax.legend(title="reported group", loc="upper left", fontsize=8.5, title_fontsize=8.5)
    vp.title(fig, f"{label}: saturated conductivity",
             f"Predicted median with its 5–95 % range. Line = perfect; bands = within\n"
             f"2× (darker) and 10× (lighter). Covariates: {ks_arm}.")
    save(fig, out, "fig4_ks")


def fig_arms(stats, out, label):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    fig.subplots_adjust(top=0.72, bottom=0.12, left=0.2, right=0.98, wspace=0.55)
    names = stats.covariates.tolist()
    yy = np.arange(len(names))[::-1]
    a = axes[0]
    a.barh(yy + 0.18, stats.exact_pct, height=0.34, color=vp.S1, label="exact class")
    a.barh(yy - 0.18, stats.group_pct, height=0.34, color=vp.S2, label="texture group")
    for yv, v1, v2 in zip(yy, stats.exact_pct, stats.group_pct):
        a.text(v1 + 1, yv + 0.18, f"{v1:.0f}", va="center", fontsize=8)
        a.text(v2 + 1, yv - 0.18, f"{v2:.0f}", va="center", fontsize=8)
    a.set_xlim(0, 100); a.set_xlabel("% of soils")
    a.legend(fontsize=8, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2)
    a.set_yticks(yy, names)
    b = axes[1]
    b.barh(yy, stats.fractions_mae, height=0.5, color=vp.S1)
    for yv, v in zip(yy, stats.fractions_mae):
        b.text(v + 0.1, yv, f"{v:.1f}", va="center", fontsize=8)
    b.set_xlabel("fractions, mean |error| (points)"); b.set_yticks(yy, [])
    c = axes[2]
    if stats.ks_n.max() > 0:
        c.barh(yy, stats.ks_median_factor, height=0.5, color=vp.S1)
        for yv, v, w, n, nt in zip(yy, stats.ks_median_factor, stats.ks_within_2x_pct,
                                   stats.ks_n, stats.n):
            extra = f"; n={n}" if n < nt else ""
            c.text(v + 0.03, yv, f"×{v:.2f}  ({w:.0f} % within 2×{extra})", va="center", fontsize=8)
        c.set_xlim(1, max(stats.ks_median_factor) * 1.6)
    else:
        c.text(0.5, 0.5, "no measured Ks in this set", ha="center", va="center",
               transform=c.transAxes, fontsize=9, color=vp.INK2)
        c.set_xticks([])
    c.set_xlabel("Ks, typical error factor"); c.set_yticks(yy, [])
    vp.title(fig, f"{label}: which covariates help",
             "The same soils predicted with each option; the best for each output is used "
             "in the other figures and the table.")
    save(fig, out, "fig5_covariate_options")


def fmt_md(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(
            f"{v:.2f}" if isinstance(v, float) else str(v) for v in r) + " |")
    return "\n".join(lines)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    name = args[0] if args else "babaeian_az"
    n_mc = int(args[1]) if len(args) > 1 else 300
    label = LABELS.get(name, name)
    out = os.path.join("figures", "external", name)
    os.makedirs(out, exist_ok=True)

    cache = os.path.join(out, "predictions_all_options.csv")
    if "--replot" in sys.argv and os.path.exists(cache):
        tg = st.load_reference_df(name).reset_index(drop=True)
        p = pd.read_csv(cache)
        p["top2"] = p.top2.map(ast.literal_eval)
    else:
        tg, p = predict(name, n_mc)
    ks_label = "measured Ks"
    if name in KS_TRUTH:
        col, ks_label = KS_TRUTH[name]
        tg["ksat_cmh"] = tg[col]
    stats = arm_stats(tg, p)
    cls, frac, ks = best_arms(stats)
    soils = soil_table(tg, p, cls, frac, ks)
    conf, gsum = group_table(soils)
    group_rows = soils[["sample", "reported class", "predicted class",
                        "reported group", "predicted group", "match"]]

    p.to_csv(cache, index=False)
    soils.to_csv(os.path.join(out, "soils.csv"), index=False)
    group_rows.to_csv(os.path.join(out, "groups.csv"), index=False)
    conf.to_csv(os.path.join(out, "groups_confusion.csv"))
    gsum.to_csv(os.path.join(out, "groups_summary.csv"), index=False)
    stats.round(3).to_csv(os.path.join(out, "stats.csv"), index=False)

    truth = tg.texture_class
    top = truth.value_counts()
    base_exact = top.iloc[0] / len(tg) * 100
    base_group = np.mean(truth.map(GROUP) == GROUP[top.index[0]]) * 100
    st_short = stats[["covariates", "exact_pct", "group_pct", "top2_pct",
                      "sand_mae", "silt_mae", "clay_mae", "fractions_mae", "fractions_n",
                      "ks_n", "ks_median_factor", "ks_bias_log10", "ks_rmse_log10",
                      "ks_within_2x_pct", "ks_within_10x_pct", "ks_cover_pct",
                      "ks_spearman"]]
    md = [f"# {label}: predicted against reported texture and Ks",
          "",
          f"{len(tg)} soils, each predicted from its measured retention curve "
          f"({CURVE_NOTE.get(name, 'points to 1500 kPa')}) with the default "
          f"reference, which does not hold "
          f"this set; n_mc={n_mc}. Always answering the most common class "
          f"('{top.index[0]}') scores {base_exact:.1f} % exact and {base_group:.1f} % "
          f"on the texture group.",
          "",
          f"Covariates used: class -- **{cls}**; fractions -- **{frac}**; "
          f"Ks -- **{ks}**. Each is the option that scored best on these same "
          f"soils (see the option table), so the figures are an upper bound. "
          f"Ks is compared with: {ks_label}.",
          "", "## Soils (sorted by clay)", "", fmt_md(soils), "",
          "## Texture groups", "", fmt_md(group_rows), "",
          "Reported (rows) against predicted (columns) group:", "",
          fmt_md(conf.reset_index()), "", fmt_md(gsum), "",
          "## Set statistics by covariate option", "",
          "Accuracy in %; fractions mean |error| in points; Ks bias and RMSE in "
          "log10 cm/h, within-factor and p5-p95 cover in %.", "",
          fmt_md(st_short), ""]
    with open(os.path.join(out, "report.md"), "w") as fh:
        fh.write("\n".join(md))
    print(f"  wrote {out}/report.md and tables")
    print(fmt_md(st_short))
    print(f"best: class={cls}, fractions={frac}, Ks={ks}")

    vp.style()
    pc = p[p.arm == cls].set_index("layer_id")
    pf = p[p.arm == frac].set_index("layer_id")
    pk = p[p.arm == ks].set_index("layer_id")
    fig_triangle(tg, pf, soils, out, frac, cls, label)
    fig_confusion(soils, conf, out, cls, label)
    fig_fractions(tg, pf, out, frac, label)
    if tg.ksat_cmh.notna().any():
        fig_ks(tg, pk, out, ks, label, ks_label)
    else:
        print("  no measured Ks in this set: Ks figure skipped")
    fig_arms(stats, out, label)


if __name__ == "__main__":
    main()
