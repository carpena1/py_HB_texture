"""Figures for "is the tool predictive for a new sensor network?".

Drawn from the caches of verify_predictive.py and verify_network.py (run
those first; nothing is recomputed here except the bootstrap intervals), in
the style of validation_plots.py. All aggregates: no individual restricted
sample is shown.

  p1_skill            texture and Ks accuracy, new and established network,
                      against depth + bulk density alone, chance and the
                      Ks replicate floor, with 95 % intervals over studies
  p2_ks_accuracy      share of soils within each factor of the measured Ks
  p3_verified_sites   gain with a network's own verified sites, unweighted
                      and weighted (verify_network.py b, d)
  p4_confidence       accuracy by the tool's own confidence signals
  p5_proximity        accuracy by distance to the nearest other-laboratory soil

Usage:  python figures/make_predictive_figures.py
"""

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

import matplotlib                                              # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                # noqa: E402
from scipy.stats import norm                                   # noqa: E402

import validation_plots as vp                                  # noqa: E402
import verify_network as vn                                    # noqa: E402
import verify_predictive as vpred                              # noqa: E402
from verify_common import GROUP                                # noqa: E402

vp.OUT = os.path.join("figures", "predictive")
NEW, EST, PHYS, BASE = vp.S1, vp.S3, vp.S2, vp.MUTED


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


def p1_skill(fr, floor, boots):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.6), gridspec_kw={"width_ratios": [3, 2.3]})
    vp.title(fig, "The curve predicts texture and Ks for a network the reference has never seen",
             "Bars: the tool with depth and bulk density, 95 % intervals over studies. "
             "Grey diamonds: depth and bulk density alone. Dotted: chance; dashed: best possible for Ks.")
    m = {k: vpred.metrics(v) for k, v in fr.items()}
    # texture
    cats = [("exact", "Exact class\n(12 classes)", "base_exact", 1 / 12),
            ("top2", "True class\nin top two", None, 2 / 12),
            ("group", "Texture group\n(4 groups)", "base_group", None)]
    w = 0.36
    for i, (key, lab, base, chance) in enumerate(cats):
        for j, (dsg, col, name) in enumerate((("source", NEW, "new network"),
                                              ("profile", EST, "established network"))):
            x = i + (j - 0.5) * w
            v = m[dsg][key] * 100
            lo, hi = np.percentile(boots[dsg][key], [2.5, 97.5]) * 100
            a1.bar(x, v, w * 0.92, color=col, label=name if i == 0 else None)
            a1.errorbar(x, v, yerr=[[v - lo], [hi - v]], color=vp.INK, lw=1, capsize=3)
            a1.text(x, hi + 1.5, f"{v:.0f}", ha="center", fontsize=8.5, color=vp.INK)
            if base:
                a1.plot(x, m[dsg][base] * 100, "D", ms=6, color=BASE, mec="white",
                        label="depth + bulk density alone" if (i == 0 and j == 0) else None)
        if chance:
            a1.plot([i - 0.42, i + 0.42], [chance * 100] * 2, ":", color=vp.INK2, lw=1.2,
                    label="chance" if i == 0 else None)
    a1.set_xticks(range(len(cats)), [c[1] for c in cats])
    a1.set_ylim(0, 100)
    a1.set_ylabel("soils predicted correctly (%)")
    a1.set_title("Texture", loc="left", fontsize=10.5, color=vp.INK)
    a1.yaxis.grid(True, color=vp.GRID, lw=0.6)
    a1.set_axisbelow(True)
    a1.legend(loc="upper left", fontsize=8.3, ncol=2)
    # Ks
    sig = floor["single"] / 0.674
    best = {"w2": floor["best_w2"], "w10": floor["best_w10"]}
    arms = [("source", "ks", NEW, "new network, neighbours"),
            ("source", "phys", PHYS, "new network, physical"),
            ("profile", "ks", EST, "established, neighbours")]
    w = 0.26
    for i, (key, lab) in enumerate((("w2", "within a\nfactor of 2"), ("w10", "within a\nfactor of 10"))):
        for j, (dsg, arm, col, name) in enumerate(arms):
            x = i + (j - 1) * w
            v = m[dsg][f"{arm}_{key}"] * 100
            lo, hi = np.percentile(boots[dsg][f"{arm}_{key}"].dropna(), [2.5, 97.5]) * 100
            a2.bar(x, v, w * 0.92, color=col, label=name if i == 0 else None)
            a2.errorbar(x, v, yerr=[[v - lo], [hi - v]], color=vp.INK, lw=1, capsize=3)
            a2.text(x, hi + 1.5, f"{v:.0f}", ha="center", fontsize=8.5, color=vp.INK)
        a2.plot(i - w, m["source"][f"base_ks_{key}"] * 100, "D", ms=6, color=BASE, mec="white")
        a2.plot(i + w, m["profile"][f"base_ks_{key}"] * 100, "D", ms=6, color=BASE, mec="white")
        a2.plot([i - 0.42, i + 0.42], [best[key] * 100] * 2, "--", color=vp.INK2, lw=1.2,
                label=f"best possible (replicate cores, ×{10 ** floor['single']:.1f})" if i == 0 else None)
    a2.set_xticks([0, 1], ["within a\nfactor of 2", "within a\nfactor of 10"])
    a2.set_ylim(0, 100)
    a2.set_ylabel("measured cores (%)")
    a2.set_title("Saturated conductivity", loc="left", fontsize=10.5, color=vp.INK)
    a2.yaxis.grid(True, color=vp.GRID, lw=0.6)
    a2.set_axisbelow(True)
    a2.legend(loc="upper left", fontsize=8.3, frameon=True, facecolor=vp.SURFACE,
               edgecolor="none", framealpha=1)
    n_s, n_k = len(fr["source"]), int((fr["source"].t_ks > 0).sum())
    fig.text(0.01, 0.01, f"{n_s} soils ({n_k} with a measured Ks) from {fr['source'].source.nunique()} "
             f"studies; new network = its study hidden, established = only its profile hidden "
             f"(verify_predictive.py).", fontsize=8, color=vp.MUTED)
    fig.tight_layout(rect=(0, 0.04, 1, 0.86))
    vp.save(fig, "p1_skill")


def p2_ks_accuracy(fr, floor):
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    vp.title(fig, "How close Ks gets: share of soils within each factor of the measured value",
             "New network unless marked. The dashed curve is the best any method can do, given how "
             "much replicate cores of one soil differ.")
    f = np.logspace(0, 2, 200)
    lf = np.log10(f)

    def ecdf(d, col):
        k = d[d.t_ks > 0]
        e = np.abs(np.log10(k[col].to_numpy(float)) - np.log10(k.t_ks.to_numpy(float)))
        return np.array([np.mean(e <= x) for x in lf]) * 100

    sig = floor["single"] / 0.674
    ax.plot(f, (2 * norm.cdf(lf / sig) - 1) * 100, "--", color=vp.INK2, lw=1.5,
            label="best possible (replicate cores)")
    for d, col, c, ls, lab in ((fr["source"], "ks_phys", PHYS, "-", "physical estimate"),
                               (fr["source"], "ks", NEW, "-", "neighbour estimate (the tool's Ks)"),
                               (fr["profile"], "ks", EST, "-.", "neighbour estimate, established network"),
                               (fr["source"], "base_ks", BASE, "-", "depth + bulk density alone"),
                               (fr["source"], "ks_median", vp.AXIS, ":", "reference median (no information)")):
        ax.plot(f, ecdf(d, col), ls, color=c, lw=2 if c in (NEW, PHYS, EST) else 1.4, label=lab)
    for x, lab in ((2, "×2"), (10, "×10")):
        ax.axvline(x, color=vp.GRID, lw=1, zorder=0)
        ax.text(x, 101, lab, ha="center", va="bottom", fontsize=9, color=vp.INK2)
    ax.set_xscale("log")
    ax.set_xlim(1, 100)
    ax.set_ylim(0, 100)
    ax.set_xticks([1, 2, 5, 10, 20, 50, 100], ["1", "2", "5", "10", "20", "50", "100"])
    ax.set_xlabel("predicted Ks within a factor of …")
    ax.set_ylabel("measured cores (%)")
    ax.yaxis.grid(True, color=vp.GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right", fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    vp.save(fig, "p2_ks_accuracy")


def p3_verified_sites(g):
    plan = vn.make_plan(g)
    est = pd.read_csv(os.path.join(vn.CACHE, "network_weighting.csv"))
    res = pd.read_csv(os.path.join(vn.CACHE, "network_pool_resid.csv"))
    per = vn.weighting_per_source(g, plan, est, res)
    b0 = per[per.k == 0].set_index("source")
    learn = pd.read_csv(os.path.join(vn.CACHE, "network_learning.csv"))
    rb = np.random.default_rng(4)
    ks_ = [k for k in vn.K_GRID if k > 0]

    def curve(arm, c):
        mean, lo, hi = [], [], []
        for k in ks_:
            x = per[(per.k == k) & (per.arm == arm)].set_index("source")
            v = (x[c] - b0.loc[x.index, c]).to_numpy()
            v = v[np.isfinite(v)] * 100
            bs = [np.mean(rb.choice(v, len(v))) for _ in range(2000)]
            mean.append(np.mean(v)); lo.append(np.percentile(bs, 2.5)); hi.append(np.percentile(bs, 97.5))
        return np.array(mean), np.array(lo), np.array(hi)

    # all of a source's other profiles (verify_network.py b)
    truth = g.set_index("layer_id")[["texture_class", "ksat_cmh"]]
    lr = learn.join(truth, on="layer_id")
    lr["ex"] = lr.cls == lr.texture_class
    lr["e"] = np.where(lr.ksat_cmh > 0, np.abs(np.log10(lr.ks / lr.ksat_cmh)), np.nan)
    pa = lr.groupby(["source", "k"]).agg(ex=("ex", "mean"),
                                         w10=("e", lambda v: np.mean(v.dropna() <= 1) if v.notna().sum() else np.nan))
    pa = pa.reset_index()
    pa["k"] = pa.k.astype(str)
    a0, aa = pa[pa.k == "0"].set_index("source"), pa[pa.k == "all"].set_index("source")
    all_ex = np.nanmean(aa.ex - a0.loc[aa.index, "ex"]) * 100
    all_w10 = np.nanmean(aa.w10 - a0.loc[aa.index, "w10"]) * 100

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.4))
    vp.title(fig, "Making a network's own verified sites count",
             "Gain over the new-network prediction as a network adds lab-verified sites (cores analysed "
             "for texture and Ks); 95 % intervals over 25 laboratories.")
    for ax, c, arms, ref, ylab, head in (
            (a1, "ex", [("base", BASE, "unweighted (the tool today)"),
                        (f"gbm x{vn.DUP}", NEW, f"own sites weighted ×{vn.DUP} in the classifier")],
             all_ex, "gain in exact class (points)", "Texture"),
            (a2, "w10", [("base", BASE, "unweighted (the tool today)"),
                         ("Ks offset", vp.S4, "local Ks offset from the verified cores"),
                         ("knn x3", PHYS, "own sites preferred by the neighbours (λ = 3)")],
             all_w10, "gain in Ks within ×10 (points)", "Saturated conductivity")):
        for arm, col, lab in arms:
            mu, lo, hi = curve(arm, c)
            ax.fill_between(ks_, lo, hi, color=col, alpha=0.15, lw=0)
            ax.plot(ks_, mu, "o-", color=col, lw=2, ms=4, label=lab)
        ax.axhline(0, color=vp.AXIS, lw=0.8)
        ax.axhline(ref, color=vp.INK2, ls=":", lw=1.2,
                   label=f"all of a laboratory's profiles, unweighted (median 111): +{ref:.0f}")
        ax.set_xscale("log")
        ax.set_xticks(ks_, [str(k) for k in ks_])
        ax.set_xlabel("verified sites of the network in the reference")
        ax.set_ylabel(ylab)
        ax.set_title(head, loc="left", fontsize=10.5, color=vp.INK)
        ax.yaxis.grid(True, color=vp.GRID, lw=0.6)
        ax.set_axisbelow(True)
        ax.legend(loc="upper left", fontsize=8.3)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    vp.save(fig, "p3_verified_sites")


def p4_confidence(fr):
    d = fr["source"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.3), gridspec_kw={"width_ratios": [1.5, 1]})
    vp.title(fig, "The tool's own signals say which answers to trust (new network)",
             "Accuracy within each group of soils; the share of soils in each group is printed under it.")
    bins = [(0, .3, "< 0.3"), (.3, .5, "0.3–0.5"), (.5, .7, "0.5–0.7"), (.7, 1.01, "≥ 0.7")]
    grp = lambda x: np.mean([GROUP[a] == GROUP[b] for a, b in zip(x.cls, x.truth)]) * 100
    for i, (lo, hi, lab) in enumerate(bins):
        x = d[(d.p1 >= lo) & (d.p1 < hi)]
        ex = np.mean(x.cls == x.truth) * 100
        a1.bar(i - 0.2, ex, 0.38, color=NEW, label="exact class" if i == 0 else None)
        a1.bar(i + 0.2, grp(x), 0.38, color=EST, label="texture group" if i == 0 else None)
        a1.text(i - 0.2, ex + 1.5, f"{ex:.0f}", ha="center", fontsize=8.5)
        a1.text(i + 0.2, grp(x) + 1.5, f"{grp(x):.0f}", ha="center", fontsize=8.5)
    shares = [np.mean((d.p1 >= lo) & (d.p1 < hi)) * 100 for lo, hi, _ in bins]
    a1.set_xticks(range(4), [f"{b[2]}\n({sh:.0f} % of soils)" for b, sh in zip(bins, shares)])
    a1.set_xlabel("probability the classifier gives its top class")
    a1.set_ylabel("correct (%)")
    a1.set_ylim(0, 100)
    a1.yaxis.grid(True, color=vp.GRID, lw=0.6)
    a1.set_axisbelow(True)
    a1.legend(loc="upper left", fontsize=8.5)
    agree = d.cls == d.knn
    k = d[d.t_ks > 0]
    close = np.abs(np.log10(k.ks / k.ks_phys)) <= 1
    e = np.abs(np.log10(k.ks / k.t_ks))
    rows = [("class: the neighbour\nvote agrees", np.mean(d[agree].cls == d[agree].truth) * 100, agree.mean(), NEW),
            ("class: it\ndisagrees", np.mean(d[~agree].cls == d[~agree].truth) * 100, (~agree).mean(), NEW),
            ("Ks: physical within\n×10 of neighbours", np.mean(e[close] <= 1) * 100, close.mean(), PHYS),
            ("Ks: further\napart", np.mean(e[~close] <= 1) * 100, (~close).mean(), PHYS)]
    for i, (lab, v, share, col) in enumerate(rows):
        a2.bar(i, v, 0.6, color=col, alpha=1 if i % 2 == 0 else 0.45)
        a2.text(i, v + 1.5, f"{v:.0f}", ha="center", fontsize=8.5)
    a2.set_xticks(range(4), [f"{r[0]}\n({r[2] * 100:.0f} % of soils)" for r in rows], fontsize=8.2)
    a2.set_ylabel("exact class / Ks within ×10 (%)")
    a2.set_ylim(0, 100)
    a2.yaxis.grid(True, color=vp.GRID, lw=0.6)
    a2.set_axisbelow(True)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    vp.save(fig, "p4_confidence")


def p5_proximity(g):
    tg = vn.targets(g)
    est = pd.read_csv(os.path.join(vn.CACHE, "network_hidden.csv"), index_col=0)
    d, k = vn.proximity_frame(g, tg, est)
    fig, ax = plt.subplots(figsize=(8, 4.4))
    vp.title(fig, "Do nearby soils from other laboratories help? (laboratory hidden)",
             "Ks tends to improve with closeness, texture does not; neither trend is significant once "
             "the laboratory is accounted for (verify_network.py a).")
    labs = vn.DIST_LABELS
    for off, col, lab, frame, col_ok in ((-0.2, NEW, "exact class", d, "ex"),
                                         (0.2, PHYS, "Ks within ×10", k, None)):
        for i, b in enumerate(labs):
            x = frame[frame.bin == b]
            if len(x) < 5:
                continue
            ok = (x[col_ok] if col_ok else (x.e <= 1)).to_numpy()
            v = ok.mean() * 100
            lo, hi = wilson(ok.sum(), len(ok))
            ax.bar(i + off, v, 0.38, color=col, label=lab if i == 1 else None)
            ax.errorbar(i + off, v, yerr=[[v - lo * 100], [hi * 100 - v]], color=vp.INK, lw=1, capsize=3)
    nd = [int((d.bin == b).sum()) for b in labs]
    nk = [int((k.bin == b).sum()) for b in labs]
    ax.set_xticks(range(len(labs)), [f"{b}\nn = {a} / {c}" for b, a, c in zip(labs, nd, nk)])
    ax.set_xlabel("distance to the nearest reference soil from another laboratory "
                  "(n = soils / soils with Ks)")
    ax.set_ylabel("correct (%)")
    ax.set_ylim(0, 100)
    ax.yaxis.grid(True, color=vp.GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    vp.save(fig, "p5_proximity")


def main():
    vp.style()
    fr, floor = vpred.frames(reuse=True)
    boots = {k: vpred.boot(v) for k, v in fr.items()}
    p1_skill(fr, floor, boots)
    p2_ks_accuracy(fr, floor)
    g = vn.load()
    p3_verified_sites(g)
    p4_confidence(fr)
    p5_proximity(g)


if __name__ == "__main__":
    main()
