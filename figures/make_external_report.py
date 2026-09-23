"""Per-soil report for an external test set: predicted against reported texture,
fractions and Ks, as tables and figures.

Every soil is predicted from its measured retention points with the default
reference (the set itself left out), as the command line would, under eight
arms: two predictor sets crossed with four covariate options. The predictor
sets are the fitted van Genuchten parameters, which the tool ships, and water
content at the fixed heads of swcc_texture.HEADS_CM. The covariate options
are the curve only, + bulk density, + organic carbon, and both; bulk density
reaches the classifier and the neighbour search, organic carbon reaches the
neighbour search only (fractions and Ks), the classifier does not take it.

For each output -- class, fractions, Ks -- the arm that scores best over the
set is reported, and named in the table's covariates columns. That choice is
made on the same soils it is scored on, so the reported numbers are an upper
bound, and doubling the arms widens that bound; all eight are in the stats
table, which is the honest comparison.

The Ks is also compared with the physical, capillary-bundle Ks and its
variants (raw, matched, air-entry capped, Peters et al. 2023, blended with the
kNN) on the default arm: ks_variants.csv and fig6.

Outputs go to figures/external/<name>/ (git-ignored while the set is not
distributed): soils.csv, groups.csv, stats.csv, ks_variants.csv, report.md
and fig1-fig6.

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

import ks_physical                                             # noqa: E402
import ksat_metrics as km                                      # noqa: E402
import swcc_texture as st                                      # noqa: E402
import validation_plots as vp                                  # noqa: E402
from verify_common import GROUP, GROUP_ORDER                   # noqa: E402
from verify_external import DATASETS                           # noqa: E402

COVARIATES = {  # name: (bulk density, organic carbon)
    "none": (False, False),
    "bulk density": (True, False),
    "organic carbon": (False, True),
    "bulk density + organic carbon": (True, True),
}
# The two predictor sets the curve can be reduced to before matching: the
# fitted van Genuchten parameters (what the tool ships) and water content at
# the fixed heads of swcc_texture.HEADS_CM.
PREDICTORS = {"vG parameters": "vg", "fixed heads": "heads"}
ARMS = {f"{pl} | {cl}": (mode, bd, oc)
        for pl, mode in PREDICTORS.items()
        for cl, (bd, oc) in COVARIATES.items()}
# Sets whose reference Ks is empty are scored against another measured Ks,
# named here, for information only.
KS_TRUTH = {"willard": ("ksat_mpd_cmh",
                        "field permeameter Ks, dry season (not in the reference)")}
LABELS = {"babaeian_az": "Arizona soils", "willard": "Laikipia soils",
          "boorowa": "Boorowa Farm soils (CSIRO, NSW)",
          "babaeian_zanjanrood": "Zanjanrood watershed soils (Iran)",
          "nj_ssir26": "New Jersey Coastal Plain soils (SSIR 26)"}
CURVE_NOTE = {"willard": "the original laboratory values, pF 0-4.2, rescaled per sample so that pF 0 sits at 95 % of porosity",
              "boorowa": "10 cm to 15 bar: suction tables and pressure plates",
              "babaeian_zanjanrood": "0-100 cm on intact cores (hanging column), 330-15,000 cm on disturbed samples (sand box and pressure plates)",
              "nj_ssir26": "saturation and 0.02-1 bar on intact cores, 2-15 bar on crushed samples; Ks on the cores"}
MATCH_COLOR = {"exact": "#2a78d6", "same group": "#e8a33d",
               "wrong group": "#d6452a"}


def predict(name, n_mc):
    tg = st.load_reference_df(name).reset_index(drop=True)
    pts = DATASETS[name].measured_points()
    ref_df = st.load_reference_df()
    ref_df = ref_df[~ref_df.layer_id.isin(tg.layer_id)].reset_index(drop=True)
    clfs = {(m, bd): st.TextureGBM(df=ref_df, feature_mode=m,
                                   covariates=["bd"] if bd else [])
            for m in set(PREDICTORS.values()) for bd in (False, True)}
    rows = []
    for arm, (mode, use_bd, use_oc) in ARMS.items():
        ref = st.GshpReference(df=ref_df, use_bd=use_bd, use_om=use_oc,
                               feature_mode=mode)
        for _, r in tg.iterrows():
            h, th = pts[r.layer_id]
            e = st.estimate(h, th, ref=ref, clf=clfs[(mode, use_bd)], n_mc=n_mc,
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
                ks_med=k["median_cmh"], ks_p5=k["p5_cmh"], ks_p95=k["p95_cmh"],
                **{f"vg_{c}": e["vg_fit"][c]
                   for c in ("thetar", "thetas", "alpha_kpa", "n", "m")},
                vg_m_free=e["vg_fit"]["m_is_free"]))
    return tg, pd.DataFrame(rows)


# The physical Ks variants scored against the tool's kNN. Each is computed on
# the tool's own van Genuchten fit of the measured points; h_min_cm is the
# air-entry cap of ks_physical.marshall_ks: none, 10 cm roughly where
# macropores begin (radius 0.15 mm), and the tool's 50 cm. "matched" divides
# by the source-balanced factor of ks_physical.MATCHING_FACTOR, refitted here
# on the reference without this set, and for each cap on its own.
PHYS_CAPS = {"no cap": None, "cap 10 cm": 10.0,
             f"cap {ks_physical.AIR_ENTRY_CM:g} cm (the tool)": ks_physical.AIR_ENTRY_CM}
DEFAULT_ARM = "fixed heads | none"


def _marshall(g, h_min_cm):
    m = np.where(g.vg_m_free.astype(bool), g.vg_m, 1.0 - 1.0 / g.vg_n)
    return ks_physical.marshall_ks(g.vg_thetar, g.vg_thetas, g.vg_alpha_kpa,
                                   g.vg_n, m=m, h_min_cm=h_min_cm)


def matching_factors(ref_df):
    """ks_physical.MATCHING_FACTOR's definition on this reference, per cap:
    the median over the sources of each one's median physics / measured."""
    d = ref_df[ref_df.ksat_cmh > 0]
    out = {}
    for cap, h in PHYS_CAPS.items():
        e = (np.log10(ks_physical.marshall_ks(d.thetar, d.thetas, d.alpha_kpa,
                                              d.n, h_min_cm=h))
             - np.log10(d.ksat_cmh.to_numpy(float)))
        ok = np.isfinite(e)
        out[cap] = 10 ** np.median(pd.Series(e[ok]).groupby(
            d.source_db.fillna("KSSL").to_numpy()[ok]).median())
    return out


def ks_variants(tg, p, factors):
    """kNN, the physics variants and half-and-half geometric blends, all on
    the default arm's fit, scored against the set's Ks."""
    g = p[p.arm == DEFAULT_ARM].set_index("layer_id").loc[tg.layer_id]
    knn = g.ks_med.to_numpy(float)
    v = {"kNN (the tool)": knn}
    for cap, h in PHYS_CAPS.items():
        raw = _marshall(g, h)
        if cap != "cap 10 cm":
            v[f"physical, {cap}, raw"] = raw
        v[f"physical, {cap}, matched"] = raw / factors[cap]
    v["Peters et al. 2023"] = ks_physical.peters_ks(g.vg_thetar, g.vg_thetas,
                                                    g.vg_alpha_kpa)
    for cap in PHYS_CAPS:
        v[f"half kNN, half {cap}, matched"] = np.sqrt(
            knn * v[f"physical, {cap}, matched"])
    obs = tg.ksat_cmh.to_numpy(float)
    rows = []
    for name, pred in v.items():
        s = km.ksat_scores(obs, pred, np.full(len(obs), np.nan),
                           np.full(len(obs), np.nan))
        ok = np.isfinite(obs) & (obs > 0) & np.isfinite(pred) & (pred > 0)
        e = np.log10(pred[ok]) - np.log10(obs[ok])
        rows.append(dict(variant=name, n=s["n"],
                         median_factor=10 ** np.median(np.abs(e)),
                         bias_log10=s["bias_log10"], within_2x_pct=s["within_2x"] * 100,
                         within_10x_pct=s["within_10x"] * 100,
                         spearman=s["spearman"]))
    return pd.DataFrame(rows), v


def arm_stats(tg, p):
    """One row of set statistics per predictor set and covariate option."""
    out = []
    for arm, g in p.groupby("arm", sort=False):
        g = g.set_index("layer_id").loc[tg.layer_id]
        truth = tg.texture_class.to_numpy()
        pred = g.pred_class.to_numpy()
        pred_set, cov = arm.split(" | ")
        row = dict(arm=arm, predictors=pred_set, covariates=cov, n=len(tg),
                   exact_pct=np.mean(pred == truth) * 100,
                   macro_f1=float(np.mean(
                       [200 * int(np.sum((pred == c) & (truth == c)))
                        / max(int(np.sum(pred == c) + np.sum(truth == c)), 1)
                        for c in sorted(set(truth))])),
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
    then RMSE. Ties keep the earlier arm, which is the vG parameters with the
    fewer covariates."""
    s = stats.reset_index(drop=True)
    cls = s.sort_values(["exact_pct", "group_pct", "top2_pct"], ascending=False,
                        kind="stable").arm.iloc[0]
    full = s[s.fractions_n == s.n]
    frac = full.sort_values("fractions_mae", kind="stable").arm.iloc[0]
    ks = s.sort_values(["ks_median_factor", "ks_rmse_log10"],
                       kind="stable").arm.iloc[0]
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
            "predictors | covariates (class)": cls,
            **{f"{x} reported": round(r[x], 1) for x in ("sand", "silt", "clay")},
            **{f"{x} predicted [p5-p95]":
               f"{f[f'pred_{x}']:.1f} [{f[f'pred_{x}_p5']:.0f}-{f[f'pred_{x}_p95']:.0f}]"
               for x in ("sand", "silt", "clay")},
            "predictors | covariates (fractions)": frac,
            "Ks measured (cm/h)": round(r.ksat_cmh, 3),
            "Ks predicted (cm/h)": round(k.ks_med, 3),
            "Ks p5-p95 (cm/h)": f"{k.ks_p5:.3g}-{k.ks_p95:.3g}",
            "Ks predicted/measured": round(k.ks_med / r.ksat_cmh, 2),
            "predictors | covariates (Ks)": ks,
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
             f"({frac_arm}).\nColour: whether the predicted class "
             f"({cls_arm}) is exact, in the same group, or wrong.")
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
             f"Predictors and covariates: {cls_arm}.")
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
             f"Predictors and covariates: {frac_arm}.")
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
             f"2× (darker) and 10× (lighter). Predictors and covariates: {ks_arm}.")
    save(fig, out, "fig4_ks")


def fig_arms(stats, out, label):
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 5.6))
    fig.subplots_adjust(top=0.82, bottom=0.1, left=0.26, right=0.98, wspace=0.55)
    names = [a.replace(" | ", ": ") for a in stats.arm]
    yy = np.arange(len(names))[::-1]
    a = axes[0]
    bars = [("exact class", stats.exact_pct, vp.S1, 0.26),
            ("texture group", stats.group_pct, vp.S2, 0.0),
            ("macro-F1", stats.macro_f1, vp.S3, -0.26)]
    for lab, v, col, off in bars:
        a.barh(yy + off, v, height=0.25, color=col, label=lab)
        for yv, x in zip(yy, v):
            a.text(x + 1, yv + off, f"{x:.0f}", va="center", fontsize=8)
    a.set_xlim(0, 100); a.set_xlabel("% of soils (macro-F1 on the same scale)")
    a.legend(fontsize=8, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=3)
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
    vp.title(fig, f"{label}: which predictors and covariates help",
             "The same soils predicted with each predictor set (the fitted vG "
             "parameters, or water content at fixed heads) and each covariate "
             "option;\nthe best for each output is used in the other figures and "
             "the table. The choice is made on these same soils, so it is an "
             "upper bound.")
    save(fig, out, "fig5_covariate_options")


def fig_ks_variants(tg, var, values, factors, out, label, ks_label):
    fig, (a, b) = plt.subplots(1, 2, figsize=(12.6, 5.8),
                               gridspec_kw=dict(width_ratios=[1.25, 1]))
    fig.subplots_adjust(top=0.8, bottom=0.11, left=0.25, right=0.98, wspace=0.3)
    yy = np.arange(len(var))[::-1]
    cols = [vp.S1 if v.startswith("kNN") else vp.S3 if v.startswith("half")
            else vp.S2 for v in var.variant]
    a.barh(yy, var.median_factor, height=0.55, color=cols)
    for yv, r in zip(yy, var.itertuples()):
        a.text(r.median_factor * 1.04, yv,
               f"×{r.median_factor:.1f}   {r.within_10x_pct:.0f} % within 10×   "
               f"bias {r.bias_log10:+.2f}   ρ {r.spearman:.2f}",
               va="center", fontsize=7.5)
    a.set_xscale("log")
    a.set_xlim(1, var.median_factor.max() * 12)
    a.set_yticks(yy, var.variant)
    a.set_xlabel("Ks, typical error factor (median |log10 error|)")
    obs = tg.ksat_cmh.to_numpy(float)
    tool = f"cap {ks_physical.AIR_ENTRY_CM:g} cm (the tool)"
    show = [("kNN (the tool)", vp.S1, "o"),
            (f"physical, {tool}, matched", vp.S2, "s"),
            (f"half kNN, half {tool}, matched", vp.S3, "^")]
    allv = np.concatenate([obs] + [values[k] for k, _, _ in show])
    allv = allv[np.isfinite(allv) & (allv > 0)]
    lim = (10 ** np.floor(np.log10(allv.min())), 10 ** np.ceil(np.log10(allv.max())))
    xx = np.array(lim)
    b.fill_between(xx, xx / 10, xx * 10, color=vp.GRID, alpha=0.55, lw=0, zorder=0)
    b.fill_between(xx, xx / 2, xx * 2, color=vp.GRID, lw=0, zorder=0)
    b.plot(xx, xx, color=vp.INK2, lw=1, zorder=1)
    for k, col, mk in show:
        b.scatter(obs, values[k], s=16, color=col, marker=mk, alpha=0.75,
                  label=k, zorder=3, lw=0)
    b.set_xscale("log"); b.set_yscale("log")
    b.set_xlim(lim); b.set_ylim(lim); b.set_aspect("equal")
    b.set_xlabel(f"{ks_label} (cm/h)"); b.set_ylabel("predicted Ks (cm/h)")
    b.legend(loc="lower right", fontsize=8)
    vp.title(fig, f"{label}: the physical Ks against the kNN",
             "Capillary bundle (Marshall 1958) on the tool's van Genuchten fit, "
             "with no pore wider than the one that empties at the cap's suction;\n"
             "raw, or divided by the matching factor (matched; refitted on the "
             "reference without this set: "
             + ", ".join(f"{k} {v:.2f}" for k, v in factors.items())
             + "),\nPeters et al. (2023), and half-and-half geometric blends "
             "with the kNN (fixed heads, no covariates).")
    save(fig, out, "fig6_ks_physics")


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
    ref_df = st.load_reference_df()
    factors = matching_factors(ref_df[~ref_df.layer_id.isin(tg.layer_id)])
    has_ks = bool((tg.ksat_cmh > 0).any())
    if has_ks:
        var, var_values = ks_variants(tg, p, factors)
        var.round(3).to_csv(os.path.join(out, "ks_variants.csv"), index=False)
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
    one = np.full(len(tg), top.index[0], dtype=object)
    base_f1 = float(np.mean(
        [200 * int(np.sum((one == c) & (truth.to_numpy() == c)))
         / max(int(np.sum(one == c) + np.sum(truth.to_numpy() == c)), 1)
         for c in sorted(set(truth))]))
    st_short = stats[["predictors", "covariates", "exact_pct", "group_pct",
                      "macro_f1", "top2_pct",
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
          f"('{top.index[0]}') scores {base_exact:.1f} % exact, {base_group:.1f} % "
          f"on the texture group and {base_f1:.1f} macro-F1 (F1 averaged over "
          f"the classes this set contains, as the overall validation reports "
          f"it).",
          "",
          f"Predictors and covariates used: class -- **{cls}**; fractions -- "
          f"**{frac}**; Ks -- **{ks}** (predictor set | covariates). Each is "
          f"the arm that scored best on these same soils (see the arm table), "
          f"so the figures are an upper bound. "
          f"Ks is compared with: {ks_label}.",
          "", "## Soils (sorted by clay)", "", fmt_md(soils), "",
          "## Texture groups", "", fmt_md(group_rows), "",
          "Reported (rows) against predicted (columns) group:", "",
          fmt_md(conf.reset_index()), "", fmt_md(gsum), "",
          "## Set statistics by predictor set and covariate option", "",
          "Accuracy in %; fractions mean |error| in points; Ks bias and RMSE in "
          "log10 cm/h, within-factor and p5-p95 cover in %.", "",
          fmt_md(st_short), ""]
    if has_ks:
        md += ["## Physical Ks against the kNN", "",
               f"On the default arm's van Genuchten fit ({DEFAULT_ARM}); see "
               f"fig6. Matching factors refitted on the reference without this "
               f"set: " + ", ".join(f"{k} {v:.2f}" for k, v in factors.items())
               + ". median_factor is the typical error factor; bias is the "
               "mean log10 error.", "", fmt_md(var), ""]
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
    if has_ks:
        print(fmt_md(var))
        fig_ks_variants(tg, var, var_values, factors, out, label, ks_label)


if __name__ == "__main__":
    main()
