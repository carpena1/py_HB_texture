"""Does the tool do worse on volcanic soils, and does telling the classifier
which soils are volcanic help?

Targets are every layer of the default reference flagged volcanic (known or
probable, prepare_volcanic.py), the andic ones reported separately, plus a
control draw of other soils (up to n_per_class per class). The Canary set,
which is not in the reference, is predicted as a new source.

Each design trains the classifier on the reference minus the held-out folds:
  profile  5 folds grouped by profile: the target's source stays in the
           reference (what a user whose data are in it would see)
  source   5 folds grouped by data source: the target's whole source is out

Arms, all from the same fitted parameters (no Monte Carlo, so the arms differ
only in their features):
  base         the four vG parameters + sample type (the shipped classifier)
  + volcanic   + volcanic parent material (known/probable 1, possible 0.5,
               unlikely 0, unknown missing)
  + andic      + andic properties (yes/likely 1, no 0, unknown missing)
  + both
A user would state the flag, so the target carries its own value.

Usage:  python verify_volcanic.py [n_per_class] [seed]

seed sets the fold split and the control draw (default 0).
"""

import sys

import numpy as np
import pandas as pd

import swcc_texture as st
from verify_common import GROUP, mcnemar

N_FOLDS = 5
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 0
VOLC = {"known": 1.0, "probable": 1.0, "possible": 0.5, "unlikely": 0.0}
ANDIC = {"yes": 1.0, "likely": 1.0, "no": 0.0}
ARMS = {"base": (), "+ volcanic": ("volcanic_code",),
        "+ andic": ("andic_code",), "+ both": ("volcanic_code", "andic_code")}


def codes(df):
    df = df.copy()
    df["volcanic_code"] = df.volcanic.map(VOLC)
    df["andic_code"] = df.andic.map(ANDIC)
    return df


def predict(clf, df):
    cols = [df[c].to_numpy(float) for c in clf.covariates]
    stype = df.sample_type.map(st.SAMPLE_TYPE_CODE).to_numpy(float)
    x = clf._features(df.thetar.to_numpy(), df.thetas.to_numpy(),
                      df.alpha_kpa.to_numpy(), df.n.to_numpy(),
                      stype=stype, covs=cols)
    return clf.clf.predict(x)


def run(ref, tg, key, extra=None):
    """Predict tg with folds of tg[key] held out of the training data."""
    out = {a: np.empty(len(tg), object) for a in ARMS}
    groups = np.array(sorted(tg[key].unique()))
    np.random.default_rng(SEED).shuffle(groups)
    for chunk in np.array_split(groups, N_FOLDS):
        test = tg[key].isin(chunk).to_numpy()
        train = ref[~ref[key].isin(chunk)]
        for arm, cov in ARMS.items():
            clf = st.TextureGBM(df=train, covariates=cov)
            out[arm][test] = predict(clf, tg[test])
    if extra is not None:              # a set outside the reference: new source
        for arm, cov in ARMS.items():
            clf = st.TextureGBM(df=ref, covariates=cov)
            out[arm] = np.concatenate([out[arm], predict(clf, extra)])
    return out


def report(title, truth, out, mask):
    t = truth[mask]
    print(f"\n  {title} (n={mask.sum()})")
    base_ok = out["base"][mask] == t
    base_g = np.array([GROUP[a] == GROUP[b] for a, b in zip(out["base"][mask], t)])
    for arm in ARMS:
        p = out[arm][mask]
        ok = p == t
        g = np.array([GROUP[a] == GROUP[b] for a, b in zip(p, t)])
        vs = ("" if arm == "base" else
              f"  exact {(ok.mean()-base_ok.mean())*100:+.1f} pp p={mcnemar(base_ok, ok):.3f}"
              f"   group {(g.mean()-base_g.mean())*100:+.1f} pp p={mcnemar(base_g, g):.3f}")
        print(f"    {arm:<11s} exact {ok.mean()*100:5.1f} %  group {g.mean()*100:5.1f} %{vs}")


def main():
    n_per_class = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    ref = codes(st.load_reference_df())
    if "volcanic" not in ref:
        raise SystemExit("no volcanic flags: run prepare_volcanic.py first")
    ref["source_db"] = ref.source_db.fillna(ref.layer_id.astype(str).str.split("_").str[0])
    ref["profile_id"] = ref.profile_id.fillna("solo_" + ref.layer_id.astype(str))
    volc = ref.volcanic.isin(["known", "probable"])
    rest = ref[~volc & ref.volcanic.isin(["unlikely"])]
    ctrl = pd.concat([g.sample(min(len(g), n_per_class), random_state=SEED)
                      for _, g in rest.groupby("texture_class")])
    tg = pd.concat([ref[volc], ctrl]).reset_index(drop=True)
    extra = None                 # Canary as a new source, when not in the reference
    if not ref.source_db.eq("Armas_Canarias").any():
        try:
            a = codes(st.load_reference_df("armas"))
            extra = a[a.volcanic.notna()].reset_index(drop=True)
        except FileNotFoundError:
            pass
    print(f"reference {len(ref)} layers; volcanic flags: "
          f"{ref.volcanic.value_counts().to_dict()}; andic: {ref.andic.value_counts().to_dict()}")
    print(f"targets: {volc.sum()} volcanic (known/probable) from "
          f"{ref[volc].source_db.nunique()} sources, {len(ctrl)} controls"
          + (f", {len(extra)} Canary (new source)" if extra is not None else ""))
    print(f"volcanic targets by source: {ref[volc].source_db.value_counts().head(8).to_dict()}")
    print(f"volcanic targets by class: {ref[volc].texture_class.value_counts().to_dict()}")

    for key, name in (("profile_id", "profile folds (source in the reference)"),
                      ("source_db", "source folds (new source)")):
        ex = extra if key == "source_db" else None
        out = run(ref, tg, key, ex)
        allt = pd.concat([tg, ex]) if ex is not None else tg
        truth = allt.texture_class.to_numpy()
        v = allt.volcanic.isin(["known", "probable"]).to_numpy()
        a = allt.andic.isin(["yes", "likely"]).to_numpy()
        print(f"\n=== {name} ===")
        report("volcanic targets", truth, out, v)
        report("andic targets", truth, out, a)
        report("volcanic, not andic", truth, out, v & ~a)
        report("control (not volcanic)", truth, out, ~v)
        if ex is not None:
            c = np.r_[np.zeros(len(tg), bool), np.ones(len(ex), bool)]
            report("Canary (not in the reference)", truth, out, c)


if __name__ == "__main__":
    main()
