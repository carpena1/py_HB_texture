"""Do the two alternative predictors survive the shipped pipeline?

verify_ivg.py compared predictor sets with one fit per soil. The two that
earned a second look are re-run here through the tool itself -- Monte Carlo
draws around the vG fit, the hybrid classifier, neighbour voting -- so the
numbers are comparable with verify_holdout.py's:

  heads      water content at HEADS_CM (saturation, 50, 100, 330, 1000,
             5000, 15000, 100000 cm) instead of the four vG parameters, in
             both the classifier and the neighbour search
  fractions  sand and clay from gradient-boosted regressors with
             silt = 100 - sand - clay (FractionGBM), against the kNN's
             neighbour-mean fractions

Both are reported with the curve alone and with depth and bulk density, on
the same targets and folds as verify_holdout.py (150 layers per class).

Designs: profile (a new site from a source the reference holds) and source
(a laboratory the reference has never seen).

Usage:  python verify_alt_predictors.py [profile|source] [n_per_class] [n_mc]
"""

import sys

import numpy as np
import pandas as pd

import swcc_texture as st
from verify_common import GROUP, mcnemar

H = np.arange(0.0, 150.0 + 0.1, 5.0)
N_FOLDS = 5
ARMS = ("vg", "heads")


def empty(n):
    return {"cls": np.empty(n, object), "knn": np.full((n, 3), np.nan),
            "gbm": np.full((n, 3), np.nan), "ks": np.full(n, np.nan)}


def run(train, tg, sel, out, n_mc):
    """Every arm on one fold's held-out targets."""
    train = train.reset_index(drop=True)
    models = {}
    for arm in ARMS:
        fm = "vg" if arm == "vg" else "heads"
        models[arm, False] = (
            st.GshpReference(df=train, feature_mode=fm),
            st.TextureGBM(df=train, feature_mode=fm),
            st.FractionGBM(df=train, feature_mode=fm))
        models[arm, True] = (
            st.GshpReference(df=train, feature_mode=fm, use_bd=True),
            st.TextureGBM(df=train, feature_mode=fm,
                          covariates=["depth_cm", "bd"]),
            st.FractionGBM(df=train, feature_mode=fm,
                           covariates=["depth_cm", "bd"]))
    for j in sel:
        r = tg.iloc[j]
        theta = st.vg_theta(H, r.thetar, r.thetas, r.alpha_kpa, r.n)
        has_bd = bool(np.isfinite(r.bd))
        for arm in ARMS:
            for cov in (False, True):
                if cov and not (has_bd and np.isfinite(r.depth_cm)):
                    continue
                ref, clf, frac = models[arm, cov]
                res = st.estimate(H, theta, ref=ref, clf=clf, clf_frac=frac,
                                  n_mc=n_mc, sample_type=r.sample_type,
                                  depth=r.depth_cm if cov else None,
                                  bulk_density=r.bd if cov else None)
                o = out[arm, cov]
                o["cls"][j] = res["texture_class"]
                o["knn"][j] = [res["fractions"][c]
                               for c in ("sand", "silt", "clay")]
                o["gbm"][j] = [res["fractions_gbm"][c]
                               for c in ("sand", "silt", "clay")]
                o["ks"][j] = res["ksat"]["median_cmh"]


def report(out, tg, truth, tfrac):
    for cov in (False, True):
        have = np.array([x is not None for x in out["vg", cov]["cls"]])
        print(f"\n=== {'curve + depth and bulk density' if cov else 'curve alone'}"
              f" (n={have.sum()}) ===")
        base = None
        for arm in ARMS:
            cls = out[arm, cov]["cls"][have]
            ex = cls == truth[have]
            gp = np.array([GROUP[a] == GROUP[b] for a, b in zip(cls, truth[have])])
            line = (f"class  {arm:6s} exact {ex.mean() * 100:5.1f}  "
                    f"group {gp.mean() * 100:5.1f}")
            if base is not None:
                line += (f"   vs vg {(ex.mean() - base[0].mean()) * 100:+5.1f} pp "
                         f"(p={mcnemar(base[0], ex):.3f})  group "
                         f"{(gp.mean() - base[1].mean()) * 100:+5.1f} "
                         f"(p={mcnemar(base[1], gp):.3f})")
            print(line)
            if arm == "vg":
                base = (ex, gp)
        print(f"{'fractions MAE %':30s}" +
              "".join(f"{c:>12s}" for c in ("sand", "silt", "clay")))
        for arm in ARMS:
            for how in ("knn", "gbm"):
                mae = [np.nanmean(np.abs(out[arm, cov][how][have, j]
                                         - tfrac[have, j])) for j in range(3)]
                label = f"{arm} / {'neighbour mean' if how == 'knn' else 'sand+clay regressors'}"
                print(f"{label:30s}" + "".join(f"{v:12.2f}" for v in mae))
        kt = tg.ksat_cmh.to_numpy()
        for arm in ARMS:
            p = out[arm, cov]["ks"]
            ok = have & np.isfinite(kt) & np.isfinite(p)
            e = np.log10(p[ok]) - np.log10(kt[ok])
            print(f"Ks     {arm:6s} {np.sqrt(np.mean(e ** 2)):.3f} dex, "
                  f"{np.mean(np.abs(e) < 1) * 100:4.1f} % within x10 "
                  f"(n={ok.sum()})")


def main():
    design = sys.argv[1] if len(sys.argv) > 1 else "profile"
    n_per_class = int(sys.argv[2]) if len(sys.argv) > 2 else 150
    n_mc = int(sys.argv[3]) if len(sys.argv) > 3 else 40

    g = st.load_reference_df("merged").reset_index(drop=True)
    g["source_db"] = g.source_db.fillna("KSSL")
    g["profile_id"] = g.profile_id.fillna("solo_" + g.layer_id.astype(str))
    tg = pd.concat([x.sample(min(len(x), n_per_class), random_state=0)
                    for _, x in g.groupby("texture_class")]
                   ).reset_index(drop=True)
    rng = np.random.default_rng(0)
    if design == "profile":
        p = np.array(sorted(tg.profile_id.unique()))
        rng.shuffle(p)
        fold_list = [("profile_id", set(x)) for x in np.array_split(p, N_FOLDS)]
    else:
        fold_list = [("source_db", {s}) for s in sorted(tg.source_db.unique())]
    print(f"reference {len(g)} layers; targets {len(tg)}; design {design}; "
          f"n_mc={n_mc}", flush=True)

    out = {(a, c): empty(len(tg)) for a in ARMS for c in (False, True)}
    for i, (key, held) in enumerate(fold_list, 1):
        sel = np.where(tg[key].isin(held))[0]
        if not len(sel):
            continue
        run(g[~g[key].isin(held)], tg, sel, out, n_mc)
        print(f"  fold {i}/{len(fold_list)}: {len(sel)} targets", flush=True)

    report(out, tg, tg.texture_class.to_numpy(),
           tg[["sand", "silt", "clay"]].to_numpy())


if __name__ == "__main__":
    main()
