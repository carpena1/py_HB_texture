"""Test the tool on an external dataset before it joins the reference.

The dataset's table (prepare_<name>.py) is predicted from its measured
retention points, as the command line would, in three situations:

  new source          the default reference without this dataset (it may
                      already have joined): the position of a new user or
                      sensor network.
  own data, site out  the dataset added to the reference with each site
                      (profile_id) held out in turn: a new site from a source
                      the reference already holds.
  own data, soil out  only the sample itself held out: the project's
                      standard leave-one-soil-out. (For willard the two
                      samples of a field are replicates at one depth, so
                      this one leaks.)

Classifier and neighbours are rebuilt from the same reduced reference in
every fold, and each sample is scored with and without depth and bulk
density. Always answering the dataset's most common class is the baseline
to beat, not chance.

Usage:  python verify_external.py willard|babaeian_zanjanrood|babaeian_az|armas|tong|boorowa
                                  [n_mc] [--truth=resin]
"""

import os
import sys
from collections import Counter

import numpy as np
import pandas as pd

import ksat_metrics as km
import prepare_armas
import prepare_babaeian_az
import prepare_babaeian_zanjanrood
import prepare_boorowa
import prepare_tong
import prepare_willard
import swcc_texture as st
from verify_common import GROUP, mcnemar

DATASETS = {"willard": prepare_willard,
            "babaeian_zanjanrood": prepare_babaeian_zanjanrood,
            "babaeian_az": prepare_babaeian_az, "armas": prepare_armas,
            "tong": prepare_tong, "boorowa": prepare_boorowa}
N_FOLDS = 20


def empty(n):
    return {"base": np.empty(n, object), "cov": np.empty(n, object),
            "knn": np.empty(n, object), "top2": np.empty(n, object),
            "frac": np.full((n, 3), np.nan), "ks": np.full((n, 3), np.nan),
            "ks_cov": np.full((n, 3), np.nan)}


def run(ref_df, tg, idx, pts, out, n_mc):
    """Predict tg rows idx from ref_df; classifier and neighbours both see
    only ref_df."""
    ref_df = ref_df.reset_index(drop=True)
    ref = st.GshpReference(df=ref_df)
    ref_bd = st.GshpReference(df=ref_df, use_bd=True)
    base = st.TextureGBM(df=ref_df)
    cov = st.TextureGBM(df=ref_df, covariates=["depth_cm", "bd"])
    for j in idx:
        r = tg.iloc[j]
        h, th = pts[r.layer_id]
        kw = dict(n_mc=n_mc, sample_type=r.sample_type)
        a = st.estimate(h, th, ref=ref, clf=base, **kw)
        # As the command line does with --depth and --bulk-density: bulk
        # density also joins the neighbour search.
        has_bd = bool(np.isfinite(r.bd))
        b = st.estimate(h, th, ref=ref_bd if has_bd else ref, clf=cov,
                        depth=r.depth_cm,
                        bulk_density=r.bd if has_bd else None, **kw)
        kb = b["ksat"]
        out["ks_cov"][j] = (kb["median_cmh"], kb["p5_cmh"], kb["p95_cmh"])
        out["base"][j] = a["texture_class"]
        out["top2"][j] = list(a["class_probabilities"])[:2]
        out["knn"][j] = next(iter(a["knn_class_probabilities"]))
        out["cov"][j] = b["texture_class"]
        out["frac"][j] = [a["fractions"][c] for c in ("sand", "silt", "clay")]
        k = a["ksat"]
        out["ks"][j] = (k["median_cmh"], k["p5_cmh"], k["p95_cmh"])


def held_out(ref, tg, key, pts, n_mc, rng):
    """Add tg to ref and predict it with folds of tg[key] held out."""
    both = pd.concat([ref, tg], ignore_index=True)
    groups = np.array(sorted(tg[key].unique()))
    rng.shuffle(groups)
    out = empty(len(tg))
    for chunk in np.array_split(groups, min(N_FOLDS, len(groups))):
        idx = np.where(tg[key].isin(chunk))[0]
        run(both[~both[key].isin(chunk)], tg, idx, pts, out, n_mc)
    return out


def report(title, tg, out):
    truth = tg.texture_class.to_numpy()
    grp = lambda p: np.array([GROUP[a] == GROUP[b] for a, b in zip(p, truth)])
    print(f"\n=== {title} (n={len(tg)}) ===")
    for a, name in (("base", "tool (curve only)"),
                    ("cov", "+ depth & bulk density"),
                    ("knn", "kNN second opinion")):
        print(f"  {name:<24s} exact {np.mean(out[a] == truth)*100:5.1f} %   "
              f"group {grp(out[a]).mean()*100:5.1f} %")
    top2 = np.mean([t in p for t, p in zip(truth, out["top2"])])
    print(f"  true class in the tool's top 2: {top2*100:.1f} %")
    ca, cb = out["base"] == truth, out["cov"] == truth
    print(f"  depth & bulk density vs curve only: "
          f"{(cb.mean()-ca.mean())*100:+.1f} pp, McNemar p={mcnemar(ca, cb):.3f}")
    agree = out["base"] == out["knn"]
    if agree.any() and (~agree).any():
        print(f"  tool and kNN agree on {agree.mean()*100:.0f} %: right "
              f"{np.mean(ca[agree])*100:.0f} % when they agree, "
              f"{np.mean(ca[~agree])*100:.0f} % when not")
    print(f"  predicted classes: {dict(Counter(out['base']).most_common())}")
    print(f"  with depth & bulk density: {dict(Counter(out['cov']).most_common())}")
    for c, g in tg.groupby("texture_class"):
        i = g.index.to_numpy()
        print(f"    {c:<16s} n={len(i):3d}  curve only "
              f"{np.mean(out['base'][i] == c)*100:5.1f} %   + depth & BD "
              f"{np.mean(out['cov'][i] == c)*100:5.1f} %")
    meas = tg[["sand", "silt", "clay"]].to_numpy(float)
    err = np.abs(out["frac"] - meas)
    print(f"  fractions, mean |error| (points): sand {err[:, 0].mean():.1f}, "
          f"silt {err[:, 1].mean():.1f}, clay {err[:, 2].mean():.1f}; "
          f"mean clay measured {meas[:, 2].mean():.0f} %, "
          f"predicted {out['frac'][:, 2].mean():.0f} %")
    obs = tg.ks_truth.to_numpy(float)
    e = np.abs(np.log10(out["ks"][:, 0]) - np.log10(obs))
    print(f"  Ks median |log10 error| {np.nanmedian(e):.2f} "
          f"(x{10 ** np.nanmedian(e):.1f}); median measured "
          f"{np.nanmedian(obs):.2f} cm/h, predicted "
          f"{np.nanmedian(out['ks'][:, 0]):.2f}")
    km.report(f"Ks -- {title}",
              [("tool", km.ksat_scores(obs, *out["ks"].T)),
               ("+ bulk density", km.ksat_scores(obs, *out["ks_cov"].T))])
    return ca


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args or args[0] not in DATASETS:
        raise SystemExit(f"usage: verify_external.py {'|'.join(DATASETS)} [n_mc]")
    name = args[0]
    n_mc = int(args[1]) if len(args) > 1 else 300
    pts = DATASETS[name].measured_points()
    tg = st.load_reference_df(name).reset_index(drop=True)
    # --truth=resin scores against an alternative texture a dataset carries
    # as sand_resin / silt_resin / clay_resin (andic soils disperse far more
    # completely with resin than with hexametaphosphate).
    alt = next((a.split("=", 1)[1] for a in sys.argv[1:]
                if a.startswith("--truth=")), None)
    if alt:
        tg = tg[tg[f"clay_{alt}"].notna()].reset_index(drop=True)
        for c in ("sand", "silt", "clay"):
            tg[c] = tg[f"{c}_{alt}"]
        tg["texture_class"] = tg[f"texture_class_{alt}"]
        name_out = f"{name}_{alt}"
    else:
        name_out = name
    # Ks is scored against ks_truth. A set that puts no Ks into the reference
    # (prepare_willard.py) is scored against its field permeameter, for
    # information only; ksat_cmh stays empty so the own-data folds see no Ks.
    tg["ks_truth"] = tg.ksat_cmh
    if tg.ksat_cmh.isna().all() and "ksat_mpd_cmh" in tg:
        tg["ks_truth"] = tg.ksat_mpd_cmh
        print("Ks scored against the field permeameter (not in the reference)")
    ref = st.load_reference_df()          # the shipped default ...
    ref = ref[~ref.layer_id.isin(tg.layer_id)]   # ... without this dataset
    truth = tg.texture_class
    top = truth.value_counts()
    print(f"{name}: {len(tg)} samples from {tg.profile_id.nunique()} sites; "
          f"classes {top.to_dict()}")
    print(f"default reference: {len(ref)} layers; n_mc={n_mc}")
    print(f"always answering '{top.index[0]}': exact "
          f"{top.iloc[0] / len(tg)*100:.1f} %, group "
          f"{np.mean(truth.map(GROUP) == GROUP[top.index[0]])*100:.1f} %")

    rng = np.random.default_rng(0)
    outs = {"new": empty(len(tg))}
    run(ref, tg, range(len(tg)), pts, outs["new"], n_mc)
    ok_new = report("new source: default reference without this dataset",
                    tg, outs["new"])
    outs["site"] = held_out(ref, tg, "profile_id", pts, n_mc, rng)
    ok_site = report("own data, each site held out", tg, outs["site"])
    outs["soil"] = held_out(ref, tg, "layer_id", pts, n_mc, rng)
    ok_soil = report("own data, only the soil held out", tg, outs["soil"])
    for label, ok in (("site out", ok_site), ("soil out", ok_soil)):
        print(f"own data ({label}) vs new source, exact class: "
              f"{(ok.mean()-ok_new.mean())*100:+.1f} pp, "
              f"McNemar p={mcnemar(ok_new, ok):.3f}")

    # Per-sample predictions stay in the git-ignored cache.
    cols = {}
    for tag, o in outs.items():
        for a in ("base", "cov", "knn"):
            cols[f"{tag}_{a}"] = o[a]
        cols[f"{tag}_clay"] = o["frac"][:, 2]
        cols[f"{tag}_ks"] = o["ks"][:, 0]
    cache = os.path.join("figures", "cache", f"external_{name_out}.csv")
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    pd.DataFrame({"layer_id": tg.layer_id, "texture_class": truth,
                  "clay": tg.clay, "ksat_cmh": tg.ks_truth, **cols}
                 ).to_csv(cache, index=False)
    print(f"wrote {cache}")


if __name__ == "__main__":
    main()
