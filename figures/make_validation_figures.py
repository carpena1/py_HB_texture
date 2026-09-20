"""Figures for the validation results (leave-one-soil-out, the standard).

Step 1 (slow, ~45 min, cached): every target layer is predicted with only
itself removed from the reference and from the classifier's training data,
by each prediction method the project has tried:

    knn         neighbour vote alone
    hybrid      the shipped default: the classifier decides the class
    + depth & bulk density, which also join the neighbour search
    heads       water content at fixed heads instead of the four vG
                parameters (swcc_texture.HEADS_CM), curve alone and with
                the covariates

and both ways of getting the particle fractions -- the neighbour mean the
tool reports and the sand/clay regressors (FractionGBM). The same run records
the nearest neighbour of each target (itself excluded) for the Cover & Hart
ceiling.
Per-sample predictions include EU-HYDI samples, so the cache lives in
figures/cache/, which is git-ignored; only aggregate figures are written.

Step 2: the figures, from the cache.

Usage:  python figures/make_validation_figures.py [--recompute]
"""

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.chdir(ROOT)

import swcc_texture as st                      # noqa: E402
from verify_common import GROUP, ORDER         # noqa: E402

CACHE = os.path.join("figures", "cache", "validation_predictions.csv")
H = np.concatenate([[0.0], np.logspace(-1, np.log10(1500.0), 25)])
N_FOLDS, N_PER_CLASS, N_MC = 5, 100, 40


def curve(r):
    return st.vg_theta(H, r.thetar, r.thetas, r.alpha_kpa, r.n)


def compute():
    df = st.load_reference_df("merged")
    df["source_db"] = df.source_db.fillna("KSSL")
    df["profile_id"] = df.profile_id.fillna("solo_" + df.layer_id.astype(str))
    df.loc[~df.bd.between(0.1, 2.3), "bd"] = np.nan
    pool = df[df.depth_cm.notna() & df.bd.notna()]
    tg = pd.concat([g.sample(min(len(g), N_PER_CLASS), random_state=0)
                    for _, g in pool.groupby("texture_class")]).reset_index(drop=True)

    # Nearest neighbour of each target in the full reference, itself excluded,
    # in the standardised coordinates the kNN uses.
    ref_all = st.GshpReference(df=df.reset_index(drop=True))
    pos = pd.Series(np.arange(len(df)), index=df.layer_id.to_numpy())
    nn = []
    for lid in tg.layer_id:
        f = ref_all.z[pos[lid]]
        d = np.sqrt(((ref_all.z - f) ** 2).sum(axis=1))
        d[pos[lid]] = np.inf
        nn.append(ref_all.classes[int(np.argmin(d))])
    tg["nn_class"] = nn

    rng = np.random.default_rng(0)
    lay = np.array(tg.layer_id)
    rng.shuffle(lay)
    cols = ["pred_base", "pred_cov", "pred_knn", "ks_med", "ks_p5", "ks_p95",
            "ks_cov_med", "ks_cov_p5", "ks_cov_p95",
            "pred_heads", "pred_heads_cov", "ks_heads_med", "ks_heads_cov_med"]
    frac_cols = [f"{how}_{c}" for how in ("fknn", "fgbm", "fgbm_cov")
                 for c in ("sand", "silt", "clay")]
    cols += frac_cols
    for c in cols:
        tg[c] = np.nan if c.startswith(("ks", "f")) else None
    for k, held in enumerate(np.array_split(lay, N_FOLDS), 1):
        held = set(held)
        train = df[~df.layer_id.isin(held)].reset_index(drop=True)
        ref = st.GshpReference(df=train)
        ref_bd = st.GshpReference(df=train, use_bd=True)
        base = st.TextureGBM(df=train, fractions=True)
        cov = st.TextureGBM(df=train, covariates=["depth_cm", "bd"],
                            fractions=True)
        ref_h = st.GshpReference(df=train, feature_mode="heads")
        ref_h_bd = st.GshpReference(df=train, feature_mode="heads",
                                    use_bd=True)
        base_h = st.TextureGBM(df=train, feature_mode="heads")
        cov_h = st.TextureGBM(df=train, feature_mode="heads",
                              covariates=["depth_cm", "bd"])
        for j in np.where(tg.layer_id.isin(held))[0]:
            r = tg.iloc[j]
            kw = dict(n_mc=N_MC, sample_type=r.sample_type)
            a = st.estimate(H, curve(r), ref=ref, clf=base, **kw)
            # What the command line does with --depth and --bulk-density.
            b = st.estimate(H, curve(r), ref=ref_bd, clf=cov,
                            depth=r.depth_cm, bulk_density=r.bd, **kw)
            tg.loc[j, ["pred_base", "pred_cov", "pred_knn"]] = [
                a["texture_class"], b["texture_class"],
                next(iter(a["knn_class_probabilities"]))]
            ks = a["ksat"]
            tg.loc[j, ["ks_med", "ks_p5", "ks_p95"]] = [
                ks["median_cmh"], ks["p5_cmh"], ks["p95_cmh"]]
            kb = b["ksat"]
            tg.loc[j, ["ks_cov_med", "ks_cov_p5", "ks_cov_p95"]] = [
                kb["median_cmh"], kb["p5_cmh"], kb["p95_cmh"]]
            # fixed-head predictors, the same two ways
            ah = st.estimate(H, curve(r), ref=ref_h, clf=base_h, **kw)
            bh = st.estimate(H, curve(r), ref=ref_h_bd, clf=cov_h,
                             depth=r.depth_cm, bulk_density=r.bd, **kw)
            tg.loc[j, ["pred_heads", "pred_heads_cov"]] = [
                ah["texture_class"], bh["texture_class"]]
            tg.loc[j, ["ks_heads_med", "ks_heads_cov_med"]] = [
                ah["ksat"]["median_cmh"], bh["ksat"]["median_cmh"]]
            # fractions: the neighbour mean the tool reports, and the
            # regressors, with and without the covariates
            for pre, res in (("fknn", a), ("fgbm", a), ("fgbm_cov", b)):
                key = "fractions" if pre == "fknn" else "fractions_gbm"
                tg.loc[j, [f"{pre}_{c}" for c in ("sand", "silt", "clay")]] = [
                    res[key][c] for c in ("sand", "silt", "clay")]
        print(f"  fold {k}/{N_FOLDS} done", flush=True)

    keep = ["layer_id", "source_db", "sample_type", "texture_class",
            "sand", "silt", "clay", "ksat_cmh", "nn_class"] + cols
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    tg[keep].to_csv(CACHE, index=False)
    print(f"wrote {CACHE} ({len(tg)} targets)")


def main():
    if "--recompute" in sys.argv or not os.path.exists(CACHE):
        compute()
    if "--compute-only" in sys.argv:
        return
    import validation_plots                    # noqa: E402  (figures/)
    validation_plots.make_all(pd.read_csv(CACHE), GROUP, ORDER)


if __name__ == "__main__":
    main()
