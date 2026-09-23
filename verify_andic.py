"""The shipped --andic option, through the full pipeline.

Targets: every andic layer of the default reference (yes or likely,
prepare_volcanic.py), a control draw of non-volcanic soils, and the Canary
set as a new source. Each target is predicted with estimate() -- classifier,
neighbours and Monte Carlo draws as on the command line -- from curves
sampled from its fitted parameters, with the flag left out and with the
target's own value stated (andic targets "yes", controls "no").

  profile  5 folds grouped by profile: the target's source stays in
  source   leave-one-source-out: only the target's own source is held out,
           as for a new user (5 folds of sources can hold out every andic
           source at once -- KSSL and Campania together are 278 of the 341
           andic layers -- and leave the flag almost nothing to learn from)

Classifier and neighbours are rebuilt from the same reduced reference in
every fold; the classifier with the flag is trained with it as a feature.
As in the tool, a target stated "yes" is matched on
swcc_texture.ANDIC_FEATURE_MODE and every other target on the default
predictors.

Each arm also reports the ensemble of the two opinions: the class
maximising the geometric mean of the classifier's and the neighbours'
probabilities (the log pool of verify_blend.py at its selected weight, 0.5,
fixed here rather than tuned on these soils).

Usage:  python verify_andic.py [n_per_class] [n_mc]
"""

import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import swcc_texture as st
from verify_common import GROUP, mcnemar

N_FOLDS = 5
H = np.concatenate([[0.0], np.logspace(-1, np.log10(1500.0), 25)])


def points(r):
    return H, st.vg_theta(H, r.thetar, r.thetas, r.alpha_kpa, r.n)


def run(ref_df, tg, key, n_mc, extra=None):
    allt = pd.concat([tg, extra], ignore_index=True) if extra is not None else tg
    n = len(allt)
    out = {a: dict(cls=np.empty(n, object), knn=np.empty(n, object),
                   ens=np.empty(n, object), frac=np.full((n, 3), np.nan), ks=np.full(n, np.nan))
           for a in ("without", "with")}
    groups = np.array(sorted(tg[key].unique()))
    np.random.default_rng(0).shuffle(groups)
    chunks = ([[g] for g in groups] if key == "source_db"
              else np.array_split(groups, N_FOLDS))
    folds = [(tg[key].isin(c).to_numpy(), ref_df[~ref_df[key].isin(c)])
             for c in chunks]
    if extra is not None:
        folds.append((np.r_[np.zeros(len(tg), bool), np.ones(len(extra), bool)],
                      ref_df))
    for test, train in folds:
        test = np.r_[test, np.zeros(n - len(test), bool)] if len(test) < n else test
        train = train.reset_index(drop=True)
        ref = st.GshpReference(df=train)
        clf0 = st.TextureGBM(df=train)
        clf1 = st.TextureGBM(df=train, covariates=["andic_code"])
        if st.ANDIC_FEATURE_MODE == ref.feature_mode:
            ref_a, clf1_a = ref, clf1
        else:
            ref_a = st.GshpReference(df=train, feature_mode=st.ANDIC_FEATURE_MODE)
            clf1_a = st.TextureGBM(df=train, covariates=["andic_code"],
                                   feature_mode=st.ANDIC_FEATURE_MODE)
        for j in np.where(test)[0]:
            r = allt.iloc[j]
            h, th = points(r)
            flag = "yes" if r.andic in ("yes", "likely") else "no"
            routed = (ref_a, clf1_a) if flag == "yes" else (ref, clf1)
            for arm, (rf, clf), a in (("without", (ref, clf0), None),
                                      ("with", routed, flag)):
                e = st.estimate(h, th, ref=rf, clf=clf, n_mc=n_mc,
                                sample_type=r.sample_type, andic=a)
                o = out[arm]
                o["cls"][j] = e["texture_class"]
                o["knn"][j] = next(iter(e["knn_class_probabilities"]))
                pg, pk = e["class_probabilities"], e["knn_class_probabilities"]
                o["ens"][j] = max(st.USDA_CLASSES, key=lambda c: (
                    np.log(max(pg.get(c, 0.0), 1e-12))
                    + np.log(max(pk.get(c, 0.0), 1e-12))))
                o["frac"][j] = [e["fractions"][c] for c in ("sand", "silt", "clay")]
                o["ks"][j] = e["ksat"]["median_cmh"]
    return allt, out


def report(title, allt, out, mask):
    t = allt.texture_class.to_numpy()[mask]
    meas = allt[["sand", "silt", "clay"]].to_numpy(float)[mask]
    obs = allt.ksat_cmh.to_numpy(float)[mask]
    print(f"\n  {title} (n={mask.sum()}, {np.isfinite(obs).sum()} with Ks)")
    res = {}
    for arm in ("without", "with"):
        o = out[arm]
        ok = o["cls"][mask] == t
        g = np.array([GROUP[a] == GROUP[b] for a, b in zip(o["cls"][mask], t)])
        kk = o["knn"][mask] == t
        ee = o["ens"][mask] == t
        eg = np.array([GROUP[a] == GROUP[b] for a, b in zip(o["ens"][mask], t)])
        err = np.abs(o["frac"][mask] - meas).mean(axis=1)
        ke = np.abs(np.log10(o["ks"][mask]) - np.log10(obs))
        res[arm] = (ok, g, kk, err, ke, ee, eg)
        print(f"    {arm:<8s} exact {ok.mean()*100:5.1f} %  group {g.mean()*100:5.1f} %  "
              f"kNN {kk.mean()*100:5.1f} %  ensemble {ee.mean()*100:5.1f} % / "
              f"group {eg.mean()*100:5.1f} %  fractions MAE {err.mean():4.1f}"
              + (f"  Ks med |log err| {np.nanmedian(ke):.2f}" if np.isfinite(ke).sum() else ""))
    (a0, g0, k0, e0, s0, _, _), (a1, g1, k1, e1, s1, x1, y1) = res["without"], res["with"]
    line = (f"    with vs without: exact {(a1.mean()-a0.mean())*100:+.1f} pp "
            f"p={mcnemar(a0, a1):.3f}, group {(g1.mean()-g0.mean())*100:+.1f} pp "
            f"p={mcnemar(g0, g1):.3f}, kNN {(k1.mean()-k0.mean())*100:+.1f} pp")
    if np.any(e1 != e0):
        line += f", fractions {(e1-e0).mean():+.2f} pts p={wilcoxon(e1, e0).pvalue:.3f}"
    fin = np.isfinite(s0) & np.isfinite(s1)
    if fin.any():
        line += f", Ks identical: {np.allclose(s0[fin], s1[fin])}"
    print(line)
    print(f"    with flag, ensemble vs classifier: exact {(x1.mean()-a1.mean())*100:+.1f} pp "
          f"p={mcnemar(a1, x1):.3f}, group {(y1.mean()-g1.mean())*100:+.1f} pp "
          f"p={mcnemar(g1, y1):.3f}")


def main():
    n_per_class = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    n_mc = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    ref_df = st.load_reference_df()
    ref_df["source_db"] = ref_df.source_db.fillna(ref_df.layer_id.astype(str).str.split("_").str[0])
    ref_df["profile_id"] = ref_df.profile_id.fillna("solo_" + ref_df.layer_id.astype(str))
    andic = ref_df.andic.isin(["yes", "likely"])
    rest = ref_df[ref_df.volcanic == "unlikely"]
    ctrl = pd.concat([g.sample(min(len(g), n_per_class), random_state=0)
                      for _, g in rest.groupby("texture_class")])
    tg = pd.concat([ref_df[andic], ctrl]).reset_index(drop=True)
    extra = None                 # Canary as a new source, when not in the reference
    if not ref_df.source_db.eq("Armas_Canarias").any():
        try:
            extra = st.load_reference_df("armas").reset_index(drop=True)
        except FileNotFoundError:
            pass
    print(f"reference {len(ref_df)}; targets {andic.sum()} andic + {len(ctrl)} "
          f"controls" + (f" + {len(extra)} Canary" if extra is not None else "")
          + f"; n_mc={n_mc}")
    for key, name in (("profile_id", "profile folds (source in the reference)"),
                      ("source_db", "leave-one-source-out (new source)")):
        allt, out = run(ref_df, tg, key, n_mc, extra if key == "source_db" else None)
        a = allt.andic.isin(["yes", "likely"]).to_numpy()
        canary = allt.source_db.eq("Armas_Canarias").to_numpy()
        print(f"\n=== {name} ===")
        report("andic targets (reference)", allt, out, a & ~canary)
        report("controls, stated not andic", allt, out, ~a)
        if canary.any():
            report("Canary Islands", allt, out, canary)


if __name__ == "__main__":
    main()
