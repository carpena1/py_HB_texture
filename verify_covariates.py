"""Do sample depth and bulk density improve the prediction for a new soil?

A user with a sensor usually knows its depth, and often the bulk density of
a core from the same spot. Each is scored against the shipped default (hybrid,
with sample type) on identical targets and folds, in two ways:

  GBM only   the covariate is an extra feature of the texture classifier;
             the kNN (fractions, Ks, intervals) is unchanged
  GBM + kNN  it is also a fifth matching dimension of the kNN, so it moves
             the fractions and Ks as well
  shipped    what the command line does with --depth and --bulk-density:
             both in the classifier, bulk density also in the kNN

Bulk density is the reference column `bd` (oven-dry, g/cm3; for EU-HYDI
from porosity where BD itself is missing). The classifier arms use the shipped
TextureGBM(covariates=...), exactly as the command line does.

Targets carry both covariates and are stratified by class. Fold designs:
  layer     only the target layers are held out (the standard validation)
  profile   profile-grouped 5-fold (same laboratories, other sites)
  lab       each target's whole laboratory held out
Earlier tests on smaller references found depth helping on average but
costing European soils 7 points, and bulk density acting as a laboratory
marker; the lab design is the one that decides.

Usage:  python verify_covariates.py [n_per_class] [n_mc] [--designs=layer,profile,lab]
"""

import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

import swcc_texture as st
from verify_common import GROUP, europe_mask, mcnemar

H = np.concatenate([[0.0], np.logspace(-1, np.log10(1500.0), 25)])
N_FOLDS = 5


def curve(row):
    return st.vg_theta(H, row.thetar, row.thetas, row.alpha_kpa, row.n)


def load():
    df = st.load_reference_df("merged")
    df["source_db"] = df.source_db.fillna("KSSL")
    df["profile_id"] = df.profile_id.fillna("solo_" + df.layer_id.astype(str))
    df.loc[~df.bd.between(0.1, 2.3), "bd"] = np.nan
    return df


ARMS = ("base", "depth GBM", "depth GBM+kNN", "BD GBM", "BD GBM+kNN",
        "depth+BD GBM", "shipped both")


def run(train, targets, idx, out, n_mc):
    train = train.reset_index(drop=True)
    ref = st.GshpReference(df=train)
    ref_d = st.GshpReference(df=train, use_depth=True)
    ref_b = st.GshpReference(df=train, use_bd=True)
    base = st.TextureGBM(df=train)
    gbm = {"depth": st.TextureGBM(df=train, covariates=["depth_cm"]),
           "bd": st.TextureGBM(df=train, covariates=["bd"]),
           "both": st.TextureGBM(df=train, covariates=["depth_cm", "bd"])}
    for j in idx:
        r = targets.iloc[j]
        kw = dict(n_mc=n_mc, sample_type=r.sample_type,
                  bulk_density=r.bd)
        th = curve(r)
        res = {
            "base": st.estimate(H, th, ref=ref, clf=base, **kw),
            "depth GBM": st.estimate(H, th, ref=ref, clf=gbm["depth"],
                                     depth=r.depth_cm, **kw),
            "depth GBM+kNN": st.estimate(H, th, ref=ref_d, clf=gbm["depth"],
                                         depth=r.depth_cm, **kw),
            "BD GBM": st.estimate(H, th, ref=ref, clf=gbm["bd"], **kw),
            "BD GBM+kNN": st.estimate(H, th, ref=ref_b, clf=gbm["bd"], **kw),
            "depth+BD GBM": st.estimate(H, th, ref=ref, clf=gbm["both"],
                                        depth=r.depth_cm, **kw),
            "shipped both": st.estimate(H, th, ref=ref_b, clf=gbm["both"],
                                        depth=r.depth_cm, **kw),
        }
        for a, x in res.items():
            f = x["fractions"]
            out[a]["cls"][j] = x["texture_class"]
            out[a]["frac"][j] = np.mean([abs(f[c] - r[c])
                                         for c in ("sand", "silt", "clay")])
            k = x["ksat"]
            out[a]["ks"][j] = (k["median_cmh"], k["p5_cmh"], k["p95_cmh"])


def predict(df, tg, folds, n_mc):
    n = len(tg)
    out = {a: {"cls": np.empty(n, object), "frac": np.full(n, np.nan),
               "ks": np.full((n, 3), np.nan)} for a in ARMS}
    for col, held in folds:
        sel = np.where(tg[col].isin(held).to_numpy())[0]
        if len(sel):
            run(df[~df[col].isin(held)], tg, sel, out, n_mc)
            print(f"  ... {len(sel)} targets done", flush=True)
    return out


def report(title, tg, out):
    truth = tg.texture_class.to_numpy()
    obs = tg.ksat_cmh.to_numpy(float)
    grp = lambda p: np.array([GROUP[x] == GROUP[t] for x, t in zip(p, truth)])
    b_ok, b_g = out["base"]["cls"] == truth, grp(out["base"]["cls"])
    e_b = np.abs(np.log10(out["base"]["ks"][:, 0]) - np.log10(obs))
    print(f"\n=== {title} (n={len(tg)}, {int(np.isfinite(obs).sum())} with Ks) ===")
    print(f"{'arm':<14s} {'exact':>6s} {'vs base':>16s} {'group':>6s} "
          f"{'vs base':>16s} {'frac MAE':>8s} {'Ks med|err|':>11s} {'p':>7s} "
          f"{'bias':>6s} {'<2x':>4s} {'<10x':>5s} {'cover':>5s} {'rho':>5s}")
    for a in ARMS:
        ok, g = out[a]["cls"] == truth, grp(out[a]["cls"])
        med, lo, hi = out[a]["ks"].T
        e = np.abs(np.log10(med) - np.log10(obs))
        bias = np.nanmean(np.log10(med) - np.log10(obs))
        m = np.isfinite(e) & np.isfinite(e_b)
        pk = (wilcoxon(e_b[m], e[m]).pvalue
              if a != "base" and m.sum() > 10 and np.any(e[m] != e_b[m]) else np.nan)
        cov = np.nanmean(np.where(np.isfinite(obs) & np.isfinite(lo),
                                  (obs >= lo) & (obs <= hi), np.nan))
        f = np.isfinite(obs) & np.isfinite(med) & (obs > 0) & (med > 0)
        rho = spearmanr(obs[f], med[f]).statistic if f.sum() > 2 else np.nan
        dc = "" if a == "base" else f"{(ok.mean()-b_ok.mean())*100:+.1f} p={mcnemar(b_ok, ok):.3f}"
        dg = "" if a == "base" else f"{(g.mean()-b_g.mean())*100:+.1f} p={mcnemar(b_g, g):.3f}"
        print(f"{a:<14s} {ok.mean()*100:5.1f}% {dc:>16s} {g.mean()*100:5.1f}% "
              f"{dg:>16s} {np.nanmean(out[a]['frac']):7.1f}% "
              f"{np.nanmedian(e):11.2f} {pk:7.3f} {bias:+6.2f} "
              f"{np.mean(e[f] <= np.log10(2))*100:3.0f}% "
              f"{np.mean(e[f] <= 1)*100:4.0f}% {cov*100:4.0f}% {rho:5.2f}")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    n_per_class = int(args[0]) if args else 100
    n_mc = int(args[1]) if len(args) > 1 else 40

    df = load()
    print(f"reference: {len(df)} layers; depth {df.depth_cm.notna().mean()*100:.1f} %, "
          f"bulk density {df.bd.notna().mean()*100:.1f} %")
    pool = df[df.depth_cm.notna() & df.bd.notna()]
    tg = pd.concat([g.sample(min(len(g), n_per_class), random_state=0)
                    for _, g in pool.groupby("texture_class")]).reset_index(drop=True)
    print(f"targets: {len(tg)} with depth and bulk density, {n_per_class}/class, "
          f"{tg.source_db.nunique()} laboratories; n_mc={n_mc}")

    rng = np.random.default_rng(0)
    lay = np.array(tg.layer_id)
    rng.shuffle(lay)
    profs = np.array(sorted(tg.profile_id.unique()))
    rng.shuffle(profs)
    designs = {
        "layer": [("layer_id", set(x)) for x in np.array_split(lay, N_FOLDS)],
        "profile": [("profile_id", set(x)) for x in np.array_split(profs, N_FOLDS)],
        "lab": [("source_db", {s}) for s in sorted(tg.source_db.unique())],
    }
    pick = next((a.split("=", 1)[1].split(",") for a in sys.argv[1:]
                 if a.startswith("--designs=")), ["profile", "lab"])
    designs = {k: v for k, v in designs.items() if k in pick}
    eu = europe_mask(tg).to_numpy()
    und = (tg.sample_type == "undisturbed").to_numpy()
    for name, folds in designs.items():
        print(f"\n--- {name} folds ---", flush=True)
        out = predict(df, tg, folds, n_mc)
        for label, m in (("all", np.ones(len(tg), bool)), ("European", eu),
                         ("non-European", ~eu), ("undisturbed", und)):
            sub = {a: {k: v[m] for k, v in d.items()} for a, d in out.items()}
            report(f"{name} folds, {label} targets", tg[m], sub)


if __name__ == "__main__":
    main()
