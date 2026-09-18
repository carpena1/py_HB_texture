"""Figures for the validation results (leave-one-soil-out, the standard).

Step 1 (slow, ~20 min, cached): every target layer is predicted with only
itself removed from the reference and from the classifier's training data,
with and without depth + bulk density (with them, bulk density also joins the
neighbour search, as on the command line). The same run records the nearest
neighbour of each target (itself excluded) for the Cover & Hart ceiling.
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
            "ks_cov_med", "ks_cov_p5", "ks_cov_p95"]
    for c in cols:
        tg[c] = np.nan if c.startswith("ks") else None
    for k, held in enumerate(np.array_split(lay, N_FOLDS), 1):
        held = set(held)
        train = df[~df.layer_id.isin(held)].reset_index(drop=True)
        ref = st.GshpReference(df=train)
        ref_bd = st.GshpReference(df=train, use_bd=True)
        base = st.TextureGBM(df=train)
        cov = st.TextureGBM(df=train, covariates=["depth_cm", "bd"])
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
        print(f"  fold {k}/{N_FOLDS} done", flush=True)

    keep = ["layer_id", "source_db", "sample_type", "texture_class",
            "ksat_cmh", "nn_class"] + cols
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
