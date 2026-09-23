"""Is the capillary-bundle Ks worth reporting beside the kNN, and should the
two be blended?

The tool's Ks is a lookup: the 30 nearest reference layers that carry a
measured Ks. ks_physical.py answers the same question from theory alone --
the retention curve is a pore-size distribution, the capillary equation turns
each suction into a pore radius, Hagen-Poiseuille makes each pore conduct
(Marshall 1958; Hillel 1980, ch. 8 and 9). It uses no reference data, so its
errors have nothing in common with the kNN's, which is what makes it useful
as a second opinion even though it is the weaker of the two on its own.

This script scores both on the same targets, with the target's whole
contributing database hidden from the kNN -- the honest setting, and the one
the physics is indifferent to, since hiding a source cannot change a number
that never looked at the reference. It answers three questions:

  1. How do the two compare, overall, per source and per texture group?
     (The physics as the tool computes it: pores capped at
     ks_physical.AIR_ENTRY_CM, raw and divided by the matching factor.)
  2. Is their agreement a useful confidence signal? (Split the targets by
     whether the two fall within a factor of 10 and score the kNN in each
     half: if agreement means anything, the kNN is better where they agree.)
  3. Does a matching factor help? Classical capillary-bundle theory scales
     its result to a measured Ks, because the curve fixes the shape of the
     pore-size distribution and not its absolute level; Marshall's form as
     used here has no such scaling. One offset per fold, fitted on the other
     sources, is that factor with nothing measured on the target soil; it is
     defined as the tool's ks_physical.MATCHING_FACTOR is, and the script
     prints the whole-reference value to check the constant.
  4. Would blending them beat the kNN? (And, last table, do other forms of
     the physics do better: the widest pore capped at an air-entry suction,
     or Peters et al. 2023, each matched and blended the same way?) The weights are fitted on the sources
     outside the scored fold, never on the fold itself, so the answer is what
     a user would get and not what the grid could be made to show. Blending
     is geometric: log K = (1-w) log K_knn + w log K_physical.

Targets are the layers with a measured Ks, stratified by class within each
source. Curves are generated from each layer's stored vG parameters and
refitted, as in the other verification scripts.

Usage:  python verify_ks_physical.py [n_per_source] [n_mc] [--reference=NAME]
"""

import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

import ks_physical
import ksat_metrics as km
import swcc_texture as st
from verify_common import COLS, GROUP, GROUP_ORDER

MIN_SOURCE = 40          # measured-Ks layers a source needs to be scored
H = np.arange(0.0, 150.0 + 0.1, 5.0)
REF_COLS = COLS + ["source_db"]
WEIGHTS = np.arange(0.0, 1.001, 0.05)
CAP_GRID = (None, 2.0, 5.0, 10.0, 20.0, 30.0, 50.0, 100.0, 200.0)  # cm


def predict(ref, clf, targets, blocked, n_mc):
    """Predict each target with its whole source hidden from the reference."""
    ref.set_excluded(blocked)
    knn, lo, hi, phys, cls, fits = [], [], [], [], [], []
    for row in targets.itertuples():
        theta = st.vg_theta(H, row.thetar, row.thetas, row.alpha_kpa, row.n)
        r = st.estimate(H, theta, ref=ref, clf=clf, n_mc=n_mc)
        k = r["ksat"]
        knn.append(k["median_cmh"]); lo.append(k["p5_cmh"]); hi.append(k["p95_cmh"])
        phys.append(k["physical_cmh"]); cls.append(r["texture_class"])
        fits.append({f"fit_{c}": r["vg_fit"][c]
                     for c in ("thetar", "thetas", "alpha_kpa", "n")})
    return (np.array(knn, float), np.array(lo, float), np.array(hi, float),
            np.array(phys, float), np.array(cls), pd.DataFrame(fits))


def blend(knn, phys, w):
    return 10.0 ** ((1.0 - w) * np.log10(knn) + w * np.log10(phys))


def logerr(pred, obs):
    return np.abs(np.log10(pred) - np.log10(obs))


def source_offsets(with_ks, fn=ks_physical.marshall_ks):
    """Each source's median log10(physical / measured) Ks, the physics taken
    from the stored vG parameters of every reference layer with a Ks."""
    phys = fn(with_ks.thetar.to_numpy(float), with_ks.thetas.to_numpy(float),
              with_ks.alpha_kpa.to_numpy(float), with_ks.n.to_numpy(float))
    e = np.log10(phys) - np.log10(with_ks.ksat_cmh.to_numpy(float))
    ok = np.isfinite(e)
    return pd.Series(e[ok]).groupby(with_ks.source_db.to_numpy()[ok]).median()


def matching_factor(d, phys, offsets):
    """The physics is uncalibrated: the curve fixes the shape of the pore-size
    distribution but not the absolute level, which is why the classical
    formulations scale to a measured Ks (Jackson 1972). This is that matching
    factor with nothing measured on the target soil -- one per fold, the
    median over the other sources of their offsets, so that each source has
    one vote whatever its size."""
    off = np.zeros(len(d))
    for s in d.source.unique():
        off[(d.source == s).to_numpy()] = -np.median(offsets.drop(s))
    return phys * 10.0 ** off, off


def nested_blend(d, obs, knn, phys):
    """Blend weights, each fitted on the sources outside the scored fold."""
    gw, grw, chosen = np.zeros(len(d)), np.zeros(len(d)), []
    for s in d.source.unique():
        te = (d.source == s).to_numpy()
        tr = ~te
        w = best_w(obs[tr], knn[tr], phys[tr])
        gw[te] = w
        per = {}
        for gp in GROUP_ORDER:
            m = tr & (d.pred_group == gp).to_numpy()
            per[gp] = best_w(obs[m], knn[m], phys[m])
            grw[te & (d.pred_group == gp).to_numpy()] = per[gp]
        chosen.append((w, per))
    return gw, grw, chosen


def best_w(obs, knn, phys):
    """The weight with the lowest median error on this (training) set."""
    if len(obs) < 10:
        return 0.0
    err = [np.median(logerr(blend(knn, phys, w), obs)) for w in WEIGHTS]
    return float(WEIGHTS[int(np.argmin(err))])


def line(label, e, n=None):
    return (f"{label:<26s} {n if n is not None else len(e):5d} "
            f"{np.median(e):7.2f} {np.mean(e <= np.log10(2))*100:5.0f}% "
            f"{np.mean(e <= 1)*100:5.0f}%")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    ref_name = next((a.split("=", 1)[1] for a in sys.argv[1:]
                     if a.startswith("--reference=")), st.DEFAULT_REFERENCE)
    n_per_source = int(args[0]) if args else 200
    n_mc = int(args[1]) if len(args) > 1 else 20

    df = st.load_reference_df(ref_name)
    df["source_db"] = df.source_db.fillna("KSSL")
    df["profile_id"] = df.profile_id.fillna("solo_" + df.layer_id.astype(str))
    ref = st.GshpReference(df=df[REF_COLS].reset_index(drop=True))
    src_all = df.source_db.to_numpy()

    with_ks = df[df.ksat_cmh.notna() & (df.ksat_cmh > 0)]
    sizes = with_ks.source_db.value_counts()
    sources = [s for s in sizes.index if sizes[s] >= MIN_SOURCE]
    print(f"reference {len(df)} layers, {len(with_ks)} with a measured Ks; "
          f"scoring the {len(sources)} sources with >= {MIN_SOURCE} of them")
    print(f"n_mc={n_mc}, up to {n_per_source} targets per source, "
          f"each source hidden from the kNN and the classifier\n")

    keep = []
    for s in sources:
        sub = with_ks[with_ks.source_db == s]
        tg = pd.concat([g.sample(min(len(g), max(1, n_per_source // 12)),
                                 random_state=0)
                        for _, g in sub.groupby("texture_class")])
        blocked = (src_all == s)
        clf = st.TextureGBM(df=df[~blocked])
        knn, lo, hi, phys, cls, fits = predict(ref, clf, tg, blocked, n_mc)
        obs = tg.ksat_cmh.to_numpy(float)
        ok = np.isfinite(knn) & np.isfinite(phys) & (knn > 0) & (phys > 0)
        keep.append(pd.DataFrame(dict(
            source=s, obs=obs[ok], knn=knn[ok], lo=lo[ok], hi=hi[ok],
            phys=phys[ok], pred_group=[GROUP[c] for c in cls[ok]],
            group=[GROUP[c] for c in tg.texture_class.to_numpy()[ok]],
            **{c: fits[c].to_numpy()[ok] for c in fits})))
        print(f"  {s:<24s} {ok.sum():4d} targets", flush=True)
    d = pd.concat(keep, ignore_index=True)

    obs, knn, phys = d.obs.to_numpy(), d.knn.to_numpy(), d.phys.to_numpy()
    offsets = source_offsets(with_ks)
    phys_m, off = matching_factor(d, phys, offsets)
    e_knn, e_phys = logerr(knn, obs), logerr(phys, obs)
    e_pm = logerr(phys_m, obs)
    nan = np.full(len(d), np.nan)
    km.report(f"KS, {len(d)} TARGETS, SOURCE-BLOCKED", [
        ("kNN (the tool)", km.ksat_scores(obs, knn, d.lo, d.hi)),
        ("physical", km.ksat_scores(obs, phys, nan, nan)),
        ("physical, matched", km.ksat_scores(obs, phys_m, nan, nan)),
    ])
    for label, e in (("physical", e_phys), ("physical, matched", e_pm)):
        print(f"  {label} - kNN: {np.median(e) - np.median(e_knn):+.2f} dex "
              f"median error, p={wilcoxon(e, e_knn).pvalue:.1e}")
    print(f"  matching factor across folds: {10 ** -off.max():.2f}x to "
          f"{10 ** -off.min():.2f}x (the physics is divided by it); "
          f"{10 ** np.median(offsets):.2f}x over all {len(offsets)} sources "
          f"with a Ks, {ks_physical.MATCHING_FACTOR:g}x in the tool; the "
          f"sources' own run from {10 ** offsets.min():.2g}x to "
          f"{10 ** offsets.max():.2g}x")

    print(f"\nBY SOURCE (median |log10 error|)")
    print(f"{'source':<24s} {'n':>5s} {'kNN':>6s} {'physical':>9s} "
          f"{'matched':>8s} {'delta':>7s}")
    print("-" * 64)
    for s_ in d.source.unique():
        m = (d.source == s_).to_numpy()
        print(f"{s_:<24s} {m.sum():5d} {np.median(e_knn[m]):6.2f} "
              f"{np.median(e_phys[m]):9.2f} {np.median(e_pm[m]):8.2f} "
              f"{np.median(e_pm[m]) - np.median(e_knn[m]):+7.2f}")

    print(f"\nBY MEASURED TEXTURE GROUP (median |log10 error|)")
    print(f"{'group':<24s} {'n':>5s} {'kNN':>6s} {'physical':>9s} "
          f"{'matched':>8s} {'delta':>7s}")
    print("-" * 64)
    for gp in GROUP_ORDER:
        m = (d.group == gp).to_numpy()
        if not m.any():
            continue
        print(f"{gp:<24s} {m.sum():5d} {np.median(e_knn[m]):6.2f} "
              f"{np.median(e_phys[m]):9.2f} {np.median(e_pm[m]):8.2f} "
              f"{np.median(e_pm[m]) - np.median(e_knn[m]):+7.2f}")

    # Does agreement mean anything? If it does, the kNN is right more often
    # where the independent opinion backs it up. Tested on the raw physics
    # (what the tool prints) and on the matched one, since a systematic
    # offset alone could make the two disagree for no informative reason.
    print(f"\nAGREEMENT AS A CONFIDENCE SIGNAL (the two within a factor of 10)")
    print(f"{'':<26s} {'n':>5s} {'kNN med':>7s} {'<2x':>6s} {'<10x':>6s}")
    print("-" * 54)
    for label, p_ in (("physical", phys), ("matched", phys_m)):
        agree = np.abs(np.log10(p_ / knn)) <= 1.0
        print(line(f"{label}: agree", e_knn[agree]))
        print(line(f"{label}: disagree", e_knn[~agree]))
        print(f"  agreement rate {agree.mean()*100:.0f} %")

    # Blending, with every weight fitted outside the scored source.
    print(f"\nBLEND, WEIGHTS FITTED ON THE OTHER SOURCES")
    print(f"{'rule':<26s} {'n':>5s} {'med|err|':>7s} {'<2x':>6s} {'<10x':>6s}")
    print("-" * 54)
    print(line("kNN (the tool)", e_knn))
    weights = {}
    for label, p_ in (("physical", phys), ("matched", phys_m)):
        gw, grw, chosen = nested_blend(d, obs, knn, p_)
        weights[label] = (grw, chosen)
        print(line(f"{label} alone", logerr(p_, obs)))
        print(line(f"{label}, fixed 0.5", logerr(blend(knn, p_, 0.5), obs)))
        for rule, w in (("global weight", gw), ("by predicted group", grw)):
            e = logerr(blend(knn, p_, w), obs)
            print(line(f"{label}, {rule}", e))
            print(f"    vs kNN: {np.median(e) - np.median(e_knn):+.2f} dex, "
                  f"p={wilcoxon(e, e_knn).pvalue:.1e}")
        print(f"    weights across the {len(chosen)} folds (median): "
              f"global {np.median([c[0] for c in chosen]):.2f}; "
              + ", ".join(f"{gp} {np.median([c[1][gp] for c in chosen]):.2f}"
                          for gp in GROUP_ORDER))

    # The weight column is the median across folds; each fold applies its
    # own, so a group can move even where the median weight is zero.
    print(f"\nBY PREDICTED GROUP, kNN against the per-group blends")
    print(f"{'group':<14s} {'n':>5s} {'kNN':>6s} {'blend':>7s} {'med w':>6s} "
          f"{'matched':>8s} {'med w':>6s}")
    print("-" * 56)
    for gp in GROUP_ORDER:
        m = (d.pred_group == gp).to_numpy()
        if not m.any():
            continue
        cells = []
        for label, p_ in (("physical", phys), ("matched", phys_m)):
            grw, chosen = weights[label]
            cells.append((np.median(logerr(blend(knn, p_, grw), obs)[m]),
                          np.median([c[1][gp] for c in chosen])))
        print(f"{gp:<14s} {m.sum():5d} {np.median(e_knn[m]):6.2f} "
              + " ".join(f"{e:7.2f} {w:6.2f}" for e, w in cells))

    # Other forms of the physics, each with its own matching factor (same
    # definition, refitted without the scored source) and a fixed half-and-
    # half blend: no cap, 2 cm as Vogel et al. (2001), 10 cm roughly where
    # macropores begin (radius 0.15 mm), and the tool's AIR_ENTRY_CM. The
    # last is chosen on these data, so its row is not an out-of-fold score;
    # the next table is.
    forms = {}
    for cap in (None, 2.0, 10.0, ks_physical.AIR_ENTRY_CM):
        lab = "no cap" if cap is None else f"cap {cap:g} cm"
        forms[f"Marshall, {lab}"] = (
            lambda tr, ts, a, n_, cap=cap:
            ks_physical.marshall_ks(tr, ts, a, n_, h_min_cm=cap))
    forms["Peters et al. 2023"] = lambda tr, ts, a, n_: ks_physical.peters_ks(tr, ts, a)
    fit = [d[f"fit_{c}"].to_numpy(float)
           for c in ("thetar", "thetas", "alpha_kpa", "n")]
    print(f"\nPHYSICAL VARIANTS (m = 1-1/n on the refitted curve; matched = "
          f"divided by the form's own factor, refitted without the scored source)")
    print(f"{'form':<34s} {'factor':>7s} {'med|err|':>8s} {'<2x':>5s} "
          f"{'<10x':>5s} {'bias':>6s} {'rho':>5s} {'vs kNN':>7s} {'p':>8s}")
    print("-" * 94)

    def row(label, pred, factor=""):
        e = logerr(pred, obs)
        b = np.mean(np.log10(pred) - np.log10(obs))
        print(f"{label:<34s} {factor:>7s} {np.median(e):8.2f} "
              f"{np.mean(e <= np.log10(2))*100:4.0f}% {np.mean(e <= 1)*100:4.0f}% "
              f"{b:+6.2f} {spearmanr(pred, obs).statistic:5.2f} "
              f"{np.median(e) - np.median(e_knn):+7.2f} "
              f"{wilcoxon(e, e_knn).pvalue if label != 'kNN (the tool)' else np.nan:8.1e}")
    row("kNN (the tool)", knn)
    for name, fn in forms.items():
        raw = fn(*fit)
        mat, off = matching_factor(d, raw, source_offsets(with_ks, fn))
        f = f"{10 ** -np.median(off):.2f}"
        row(f"{name}, raw", raw)
        row(f"{name}, matched", mat, f)
        row(f"  half kNN, half matched", blend(knn, mat, 0.5))

    # The cap itself chosen out of fold: for each scored source, the cap
    # whose matched physics has the lowest median error on the other
    # sources' targets, each of those matched with a factor that leaves out
    # both its own source and the scored one. So neither the cap nor any
    # factor has seen the scored source.
    print(f"\nAIR-ENTRY CAP CHOSEN ON THE OTHER SOURCES (grid "
          f"{', '.join('none' if c is None else f'{c:g} cm' for c in CAP_GRID)})")
    raws, offs = {}, {}
    for c in CAP_GRID:
        fn = lambda tr, ts, a, n_, c=c: ks_physical.marshall_ks(tr, ts, a, n_,
                                                                 h_min_cm=c)
        raws[c], offs[c] = fn(*fit), source_offsets(with_ks, fn)
    src = d.source.to_numpy()
    nested, chosen = np.full(len(d), np.nan), {}
    for s_ in d.source.unique():
        te = src == s_
        err = {}
        for c in CAP_GRID:
            o = offs[c].drop(s_)
            off_tr = np.array([-np.median(o.drop(x)) for x in src[~te]])
            err[c] = np.median(logerr(raws[c][~te] * 10.0 ** off_tr, obs[~te]))
        best = min(CAP_GRID, key=lambda c: err[c])
        chosen[s_] = best
        nested[te] = raws[best][te] * 10.0 ** -np.median(offs[best].drop(s_))
    print(f"{'form':<34s} {'factor':>7s} {'med|err|':>8s} {'<2x':>5s} "
          f"{'<10x':>5s} {'bias':>6s} {'rho':>5s} {'vs kNN':>7s} {'p':>8s}")
    print("-" * 94)
    row("kNN (the tool)", knn)
    row("cap chosen per fold, matched", nested)
    row("  half kNN, half matched", blend(knn, nested, 0.5))
    counts = pd.Series(["none" if c is None else f"{c:g} cm"
                        for c in chosen.values()]).value_counts()
    print(f"  cap chosen across the {len(chosen)} folds: "
          + ", ".join(f"{k} x{v}" for k, v in counts.items()))
    for c in CAP_GRID:
        lab = "none" if c is None else f"{c:g} cm"
        print(f"  whole reference, cap {lab}: factor "
              f"{10 ** np.median(offs[c]):.2f}")


if __name__ == "__main__":
    main()
