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
Targets are the same 1,711 layers as verify_hybrid.py (150 per class).

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
    for j in sel:
        r = tg.iloc[j]
        theta = st.vg_theta(H, r.thetar, r.thetas, r.alpha_kpa, r.n)
        res = st.estimate(H, theta, ref=ref, n_mc=n_mc, clf=clf,
                          sample_type=r.sample_type)
        out["cls"][j] = res["texture_class"]
        out["knn"][j] = next(iter(res["knn_class_probabilities"]))
        k = res["ksat"]
        out["ks"][j] = (k["median_cmh"], k["p5_cmh"], k["p95_cmh"])


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
               "ks": np.full((n, 3), np.nan)}
        for col, held in folds:
            sel = np.where(tg[col].isin(held).to_numpy())[0]
            if len(sel):
                run(g[~g[col].isin(held)], tg, sel, out, n_mc)
        res[name] = out
        print(f"  {name} done", flush=True)

    grp = lambda p: np.array([GROUP[x] == GROUP[t] for x, t in zip(p, truth)])
    print(f"\n{'held out':<9s} {'exact':>6s} {'group':>6s} {'kNN':>6s} "
          f"{'Ks n':>5s} {'med|err|':>8s} {'<2x':>4s} {'<10x':>5s} {'cover':>5s} {'rho':>5s}")
    for name, o in res.items():
        med, lo, hi = o["ks"].T
        ok = np.isfinite(obs) & np.isfinite(med) & (obs > 0) & (med > 0)
        e = np.abs(np.log10(med[ok]) - np.log10(obs[ok]))
        cov = np.mean((obs[ok] >= lo[ok]) & (obs[ok] <= hi[ok]))
        rho = spearmanr(obs[ok], med[ok]).statistic
        print(f"{name:<9s} {(o['cls'] == truth).mean()*100:5.1f}% "
              f"{grp(o['cls']).mean()*100:5.1f}% {(o['knn'] == truth).mean()*100:5.1f}% "
              f"{ok.sum():5d} {np.median(e):8.2f} {np.mean(e <= np.log10(2))*100:3.0f}% "
              f"{np.mean(e <= 1)*100:4.0f}% {cov*100:4.0f}% {rho:5.2f}")
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
