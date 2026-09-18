"""How much does it matter what is held out when a test soil is hidden?

Three designs on identical targets, each removing only what it must from both
the neighbour reference and the classifier's training data:

  layer    only the target layer itself -- a new depth at a site whose other
           horizons are already in the reference
  profile  the target's whole profile -- a new site from a data source the
           reference already holds
  source   the target's whole contributing database or laboratory -- a new
           site from a source the reference has never seen

The tool itself is the shipped default (hybrid, each target's own sample type).
Targets are the same 1,711 layers as verify_hybrid.py (150 per class). Each
design is also run as the command line does with --depth and --bulk-density
(classifier trained with both, bulk density in the neighbour search), on the
targets that carry both, and compared with the curve alone on those targets.

Usage:  python verify_holdout.py [n_per_class] [n_mc]
"""

import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import swcc_texture as st
from verify_common import GROUP, mcnemar

H = np.arange(0.0, 150.0 + 0.1, 5.0)
N_FOLDS = 5


def run(train, tg, sel, out, n_mc):
    train = train.reset_index(drop=True)
    ref = st.GshpReference(df=train)
    clf = st.TextureGBM(df=train)
    ref_bd = st.GshpReference(df=train, use_bd=True)
    clf_cov = st.TextureGBM(df=train, covariates=["depth_cm", "bd"])
    for j in sel:
        r = tg.iloc[j]
        theta = st.vg_theta(H, r.thetar, r.thetas, r.alpha_kpa, r.n)
        res = st.estimate(H, theta, ref=ref, n_mc=n_mc, clf=clf,
                          sample_type=r.sample_type)
        out["cls"][j] = res["texture_class"]
        out["knn"][j] = next(iter(res["knn_class_probabilities"]))
        k = res["ksat"]
        out["ks"][j] = (k["median_cmh"], k["p5_cmh"], k["p95_cmh"])
        if np.isfinite(r.depth_cm) and np.isfinite(r.bd):
            res = st.estimate(H, theta, ref=ref_bd, n_mc=n_mc, clf=clf_cov,
                              sample_type=r.sample_type, depth=r.depth_cm,
                              bulk_density=r.bd)
            out["cls_cov"][j] = res["texture_class"]
            k = res["ksat"]
            out["ks_cov"][j] = (k["median_cmh"], k["p5_cmh"], k["p95_cmh"])


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    n_per_class = int(args[0]) if args else 150
    n_mc = int(args[1]) if len(args) > 1 else 40

    g = st.load_reference_df("merged").reset_index(drop=True)
    g["source_db"] = g.source_db.fillna("KSSL")
    g["profile_id"] = g.profile_id.fillna("solo_" + g.layer_id.astype(str))
    tg = pd.concat([x.sample(min(len(x), n_per_class), random_state=0)
                    for _, x in g.groupby("texture_class")]).reset_index(drop=True)
    print(f"reference {len(g)} layers; targets {len(tg)} from "
          f"{tg.source_db.nunique()} sources; n_mc={n_mc}\n")

    rng = np.random.default_rng(0)
    lay = np.array(tg.layer_id); rng.shuffle(lay)
    profs = np.array(sorted(tg.profile_id.unique())); rng.shuffle(profs)
    designs = {
        "layer": [("layer_id", set(x)) for x in np.array_split(lay, N_FOLDS)],
        "profile": [("profile_id", set(x)) for x in np.array_split(profs, N_FOLDS)],
        "source": [("source_db", {s}) for s in sorted(tg.source_db.unique())],
    }
    truth = tg.texture_class.to_numpy()
    obs = tg.ksat_cmh.to_numpy(float)
    res = {}
    for name, folds in designs.items():
        n = len(tg)
        out = {"cls": np.empty(n, object), "knn": np.empty(n, object),
               "ks": np.full((n, 3), np.nan), "cls_cov": np.empty(n, object),
               "ks_cov": np.full((n, 3), np.nan)}
        for col, held in folds:
            sel = np.where(tg[col].isin(held).to_numpy())[0]
            if len(sel):
                run(g[~g[col].isin(held)], tg, sel, out, n_mc)
        res[name] = out
        print(f"  {name} done", flush=True)

    grp = lambda p, m=slice(None): np.array(
        [GROUP[x] == GROUP[t] for x, t in zip(p, truth[m])])

    def ks_line(ks, m):
        med, lo, hi = ks.T
        ok = m & np.isfinite(obs) & np.isfinite(med) & (obs > 0) & (med > 0)
        e = np.abs(np.log10(med[ok]) - np.log10(obs[ok]))
        cov = np.mean((obs[ok] >= lo[ok]) & (obs[ok] <= hi[ok]))
        rho = spearmanr(obs[ok], med[ok]).statistic
        return (f"{ok.sum():5d} {np.median(e):8.2f} {np.mean(e <= np.log10(2))*100:3.0f}% "
                f"{np.mean(e <= 1)*100:4.0f}% {cov*100:4.0f}% {rho:5.2f}")

    print(f"\n{'held out':<9s} {'exact':>6s} {'group':>6s} {'kNN':>6s} "
          f"{'Ks n':>5s} {'med|err|':>8s} {'<2x':>4s} {'<10x':>5s} {'cover':>5s} {'rho':>5s}")
    every = np.ones(len(tg), bool)
    for name, o in res.items():
        print(f"{name:<9s} {(o['cls'] == truth).mean()*100:5.1f}% "
              f"{grp(o['cls']).mean()*100:5.1f}% {(o['knn'] == truth).mean()*100:5.1f}% "
              + ks_line(o["ks"], every))

    has = np.array([c is not None for c in res["layer"]["cls_cov"]])
    print(f"\nwith --depth and --bulk-density, on the {has.sum()} targets carrying both "
          f"(curve alone on the same targets above each):")
    for name, o in res.items():
        for label, c, ks in (("curve", o["cls"], o["ks"]),
                             ("+ d & BD", o["cls_cov"], o["ks_cov"])):
            print(f"{name + ' ' + label:<18s} {(c[has] == truth[has]).mean()*100:5.1f}% "
                  f"{grp(c[has], has).mean()*100:5.1f}%        " + ks_line(ks, has))
        a, b = o["cls"][has] == truth[has], o["cls_cov"][has] == truth[has]
        ga, gb = grp(o["cls"][has], has), grp(o["cls_cov"][has], has)
        print(f"  {name}: + depth & BD {(b.mean() - a.mean())*100:+.1f} pp exact "
              f"p={mcnemar(a, b):.4f}, group {(gb.mean() - ga.mean())*100:+.1f} pp "
              f"p={mcnemar(ga, gb):.4f}")
    for a, b in (("layer", "profile"), ("profile", "source")):
        ca, cb = res[a]["cls"] == truth, res[b]["cls"] == truth
        print(f"  {b} vs {a}: {(cb.mean() - ca.mean())*100:+.1f} pp exact, "
              f"McNemar p={mcnemar(ca, cb):.4f}")

    # Targets whose profile has other horizons: do those siblings drive the
    # layer-level result?
    sibs = g.groupby("profile_id").layer_id.transform("size")
    multi = tg.profile_id.map(g.groupby("profile_id").size()).to_numpy() > 1
    share = np.array([
        (g[(g.profile_id == p) & (g.layer_id != lid)].texture_class == c).any()
        for p, lid, c in zip(tg.profile_id, tg.layer_id, truth)])
    for label, m in (("targets with a same-class sibling", share),
                     ("targets with siblings, none same class", multi & ~share),
                     ("targets without siblings", ~multi)):
        if m.any():
            print(f"\n{label} (n={m.sum()}): exact " + ", ".join(
                f"{k} {(res[k]['cls'][m] == truth[m]).mean()*100:.1f} %"
                for k in res))


if __name__ == "__main__":
    main()
