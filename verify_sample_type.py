"""Does knowing how a sample was prepared help? (sample_type)

Every reference layer carries sample_type -- undisturbed, disturbed or
unknown -- set by how the WET end of its retention curve was measured (see
the prepare_*.py scripts). Near saturation an undisturbed core keeps its
macropores and a repacked one does not, so the type should matter most for
Ks. Three uses are scored against the tool without it, on identical folds:

 1. Ks from neighbours of the target's own type only (what estimate() does
    when given sample_type; class vote and fractions still use every row).
 2. sample_type as a fifth GBM feature (0/1, missing for unknown), for the
    hybrid's texture class. 1 + 2 together are the shipped default.
 3. The kNN matched on its own type for EVERYTHING, class included -- tested
    and rejected (it loses accuracy against an unseen laboratory).

Targets are layers of known type with a measured Ks, stratified by type and
texture class, so the same soils answer both the Ks and the texture question.
Two fold designs:
  profile   profile-grouped 5-fold, the in-distribution view;
  lab       each target's whole contributing laboratory (source_db) is held
            out of the reference and the GBM, the view a new user is in.
Only 45 reference layers with a measured Ks are disturbed (ETH literature,
UNSODA, AfSPDB), so disturbed targets are few and have few disturbed Ks
neighbours; the mean number of Ks-bearing neighbours is printed so that shows.

Usage:  python verify_sample_type.py [n_per_class] [n_mc] [--reference=NAME]
"""

import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import ksat_metrics as km
import swcc_texture as st
from verify_common import mcnemar
from verify_groups import GROUP

H = np.concatenate([[0.0], np.logspace(-1, np.log10(1500.0), 25)])
N_FOLDS = 5
TYPES = ("undisturbed", "disturbed")


def curve(row):
    return st.vg_theta(H, row.thetar, row.thetas, row.alpha_kpa, row.n)


ARMS = ("kNN", "kNN same type", "hybrid", "hybrid +type")


def run(train, targets, idx, out, n_mc):
    """Predict targets.iloc[idx] from `train` under every arm."""
    ref = st.GshpReference(df=train.reset_index(drop=True))
    own = {t: st.GshpReference(df=train[train.sample_type == t]
                               .reset_index(drop=True)) for t in TYPES}
    plain = st.TextureGBM(df=train.reset_index(drop=True),
                          use_sample_type=False)
    typed = st.TextureGBM(df=train.reset_index(drop=True))
    for j in idx:
        row = targets.iloc[j]
        a = st.estimate(H, curve(row), ref=ref, n_mc=n_mc, clf=plain)
        b = st.estimate(H, curve(row), ref=ref, n_mc=n_mc, clf=typed,
                        sample_type=row.sample_type)
        c = st.estimate(H, curve(row), ref=own[row.sample_type], n_mc=n_mc)
        out["cls"]["kNN"][j] = next(iter(a["knn_class_probabilities"]))
        out["cls"]["hybrid"][j] = a["texture_class"]
        out["cls"]["kNN same type"][j] = c["texture_class"]
        out["cls"]["hybrid +type"][j] = b["texture_class"]
        for arm, r in (("kNN", a), ("kNN same type", b)):
            k = r["ksat"]
            out["ks"][arm][j] = (k["median_cmh"], k["p5_cmh"], k["p95_cmh"],
                                 k["n_neighbors_with_ksat"])


def predict(ref_df, targets, folds, n_mc):
    n = len(targets)
    out = {"cls": {a: np.empty(n, dtype=object) for a in ARMS},
           "ks": {a: np.full((n, 4), np.nan) for a in ("kNN", "kNN same type")}}
    for held_col, held in folds:
        sel = np.where(targets[held_col].isin(held).to_numpy())[0]
        if len(sel):
            run(ref_df[~ref_df[held_col].isin(held)], targets, sel, out, n_mc)
            print(f"  ... fold {sorted(held)[:2]}{'...' if len(held) > 2 else ''}"
                  f" ({len(sel)} targets)", flush=True)
    return out


def report(title, tg, out):
    truth = tg.texture_class.to_numpy()
    obs = tg.ksat_cmh.to_numpy(float)
    grp = lambda p: np.array([GROUP[x] == GROUP[t] for x, t in zip(p, truth)])
    print(f"\n=== {title} (n={len(tg)}) ===")
    print(f"{'arm':<16s} {'exact':>7s} {'group':>7s}")
    for a in ARMS:
        print(f"{a:<16s} {(out['cls'][a] == truth).mean()*100:6.1f}% "
              f"{grp(out['cls'][a]).mean()*100:6.1f}%")
    for a, b in (("kNN", "kNN same type"), ("hybrid", "hybrid +type")):
        ca, cb = out["cls"][a] == truth, out["cls"][b] == truth
        ga, gb = grp(out["cls"][a]), grp(out["cls"][b])
        print(f"  {b} vs {a}: class {(cb.mean()-ca.mean())*100:+.1f}pp "
              f"p={mcnemar(ca, cb):.4f}   group {(gb.mean()-ga.mean())*100:+.1f}pp"
              f" p={mcnemar(ga, gb):.4f}")
    ks = [(a, km.ksat_scores(obs, *out["ks"][a][:, :3].T))
          for a in ("kNN", "kNN same type")]
    km.report(f"KSAT -- {title}", ks)
    for a in ("kNN", "kNN same type"):
        print(f"  {a}: mean Ks-bearing neighbours of k=30: "
              f"{np.nanmean(out['ks'][a][:, 3]):.1f}")
    ea = np.abs(np.log10(out["ks"]["kNN"][:, 0]) - np.log10(obs))
    eb = np.abs(np.log10(out["ks"]["kNN same type"][:, 0]) - np.log10(obs))
    ok = np.isfinite(ea) & np.isfinite(eb)
    if ok.sum() > 10:
        print(f"  |log10 error| same type vs all: median "
              f"{np.median(eb[ok]):.2f} vs {np.median(ea[ok]):.2f}, "
              f"Wilcoxon p={wilcoxon(ea[ok], eb[ok]).pvalue:.4f} (n={ok.sum()})")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    ref_name = next((a.split("=", 1)[1] for a in sys.argv[1:]
                     if a.startswith("--reference=")), st.DEFAULT_REFERENCE)
    n_per_class = int(args[0]) if args else 40
    n_mc = int(args[1]) if len(args) > 1 else 40

    df = st.load_reference_df(ref_name)
    df["source_db"] = df.source_db.fillna("KSSL")
    df["profile_id"] = df.profile_id.fillna("solo_" + df.layer_id.astype(str))
    df["sample_type"] = df.sample_type.fillna("unknown")
    print(f"reference {ref_name}: {len(df)} layers; sample type "
          f"{df.sample_type.value_counts().to_dict()}")
    has_ks = df.ksat_cmh.notna() & (df.ksat_cmh > 0)
    print("with measured Ks: "
          f"{df[has_ks].sample_type.value_counts().to_dict()}")

    pool = df[has_ks & df.sample_type.isin(TYPES)]
    tg = pd.concat([g.sample(min(len(g), n_per_class), random_state=0)
                    for _, g in pool.groupby(["sample_type", "texture_class"])]
                   ).reset_index(drop=True)
    print(f"targets: {len(tg)} ({tg.sample_type.value_counts().to_dict()}), "
          f"{n_per_class}/class/type; n_mc={n_mc}")
    print("target laboratories: "
          f"{tg.groupby('sample_type').source_db.nunique().to_dict()}\n")

    rng = np.random.default_rng(0)
    profs = np.array(sorted(tg.profile_id.unique()))
    rng.shuffle(profs)
    designs = {
        "profile": [("profile_id", set(x)) for x in np.array_split(profs, N_FOLDS)],
        "lab": [("source_db", {s}) for s in sorted(tg.source_db.unique())],
    }
    for name, folds in designs.items():
        print(f"--- {name} folds ---", flush=True)
        out = predict(df, tg, folds, n_mc)
        for t in TYPES + ("all",):
            m = np.ones(len(tg), bool) if t == "all" else (tg.sample_type == t).to_numpy()
            sub = {"cls": {a: v[m] for a, v in out["cls"].items()},
                   "ks": {a: v[m] for a, v in out["ks"].items()}}
            report(f"{name} folds, {t} targets", tg[m], sub)
        if name == "lab":
            print("\n  lab folds, Ks by laboratory (median |log10 error|, "
                  "all vs same type, n):")
            obs = tg.ksat_cmh.to_numpy(float)
            for s, g in tg.groupby("source_db"):
                i = g.index.to_numpy()
                ea = np.abs(np.log10(out["ks"]["kNN"][i, 0]) - np.log10(obs[i]))
                eb = np.abs(np.log10(out["ks"]["kNN same type"][i, 0])
                            - np.log10(obs[i]))
                print(f"    {s:<24s} {g.sample_type.iloc[0]:<12s} "
                      f"{np.nanmedian(ea):5.2f} {np.nanmedian(eb):5.2f}  "
                      f"{len(g):4d}")


if __name__ == "__main__":
    main()
