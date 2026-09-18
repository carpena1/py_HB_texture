"""What does EU-HYDI add?

EU-HYDI roughly quadruples the European part of the reference. Three
questions, each reported for texture (exact class, group) and Ksat, for both
the kNN and the hybrid model:

 1. Does it help the European soils ALREADY in the reference, and does it
    leave non-European soils alone? Targets come from the public tables only, so the
    gain cannot come from EU-HYDI predicting itself.

 2. How well is a European soil from a laboratory the reference has never
    seen predicted, and how much do OTHER European laboratories help it?
    EU-HYDI targets are predicted from the public tables alone, then from those plus
    EU-HYDI with the target's whole contributor held out (contributor-grouped
    folds). The second arm is the honest number for a new European user.

 3. How much of the in-source score is protocol leakage? The same targets
    with only their own profile held out, so their laboratory stays in.

Profile- or contributor-grouped folds throughout, so the kNN reference and the
GBM's training frame are always the same frame and neither sees the target.

RESTRICTED DATA: needs data/euhydi_reference.csv, built locally with
prepare_euhydi.py; it is not distributed with the repository.

Usage:  python verify_euhydi.py [n_per_class] [n_mc]
"""

import sys

import numpy as np
import pandas as pd

import swcc_texture as st
from verify_hohenbrink import N_FOLDS, block, predict
from verify_region import europe_mask


def folds_of(ids, k, seed):
    ids = np.array(sorted(set(ids)), dtype=object)
    np.random.default_rng(seed).shuffle(ids)
    return [set(x) for x in np.array_split(ids, k)]


def stratified(pool, n):
    return pd.concat([g.sample(min(len(g), n), random_state=0)
                      for _, g in pool.groupby("texture_class")])


def main():
    n_per_class = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    n_mc = int(sys.argv[2]) if len(sys.argv) > 2 else 30

    base = st.load_reference_df("public")
    eu = st.load_reference_df("euhydi")
    both = pd.concat([base, eu], ignore_index=True)
    for d in (base, eu, both):
        d["profile_id"] = d.profile_id.fillna("solo_" + d.layer_id.astype(str))

    em = europe_mask(base)
    print(f"public {len(base)}  EU-HYDI {len(eu)}  merged {len(both)}")
    print(f"European layers: {int(em.sum())} in public, "
          f"{int(europe_mask(both).sum())} with EU-HYDI")
    print(f"EU-HYDI contributors: {eu.source_db.nunique()}")
    print(f"targets {n_per_class}/class; {N_FOLDS}-fold grouped; n_mc={n_mc}\n")

    # 1 -- soils already in the reference
    for label, mask in [("EUROPEAN soils already in the public tables", em),
                        ("NON-EUROPEAN soils", ~em)]:
        tg = stratified(base[mask], n_per_class)
        folds = folds_of(tg.profile_id, N_FOLDS, 0)
        res = {}
        for tag, hyb in (("kNN", False), ("hybrid", True)):
            res[f"{tag} public"] = predict(base, tg, n_mc, hyb, folds)
            res[f"{tag} +EU-HYDI"] = predict(both, tg, n_mc, hyb, folds)
        block(label, tg.texture_class.to_numpy(), tg.ksat_cmh.to_numpy(float),
              res)

    # 2 and 3 -- EU-HYDI soils themselves
    tg = stratified(eu, n_per_class)
    truth, obs = tg.texture_class.to_numpy(), tg.ksat_cmh.to_numpy(float)
    by_src = eu.groupby("source_db").profile_id.apply(set)
    cfolds = [set().union(*(by_src[s] for s in f))
              for f in folds_of(tg.source_db, N_FOLDS, 1)]
    pfolds = folds_of(tg.profile_id, N_FOLDS, 2)

    arms = {}
    for tag, hyb in (("kNN", False), ("hybrid", True)):
        arms[f"{tag}, no EU-HYDI"] = predict(base, tg, n_mc, hyb)
        arms[f"{tag}, + other EU labs"] = predict(both, tg, n_mc, hyb, cfolds)
        arms[f"{tag}, own lab present"] = predict(both, tg, n_mc, hyb, pfolds)

    block("EU-HYDI soils from an UNSEEN laboratory", truth, obs,
          {k: arms[k] for k in ("kNN, no EU-HYDI", "kNN, + other EU labs",
                                "hybrid, no EU-HYDI", "hybrid, + other EU labs")})
    block("EU-HYDI soils: protocol leakage (own lab held out vs present)",
          truth, obs,
          {k: arms[k] for k in ("kNN, + other EU labs", "kNN, own lab present",
                                "hybrid, + other EU labs",
                                "hybrid, own lab present")})


if __name__ == "__main__":
    main()
