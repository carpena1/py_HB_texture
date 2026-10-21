"""Does KSSL's sample_type label change anything?

KSSL rows are labelled "undisturbed", but their own provenance says the wet
end is on natural clods and the dry end on material sieved to <2 mm. The
label can only act in two places, and one is ruled out by construction:

  Ks    KSSL carries no measured Ksat, so its rows are never Ks donors
        (asserted below). Nothing to test.
  class sample_type is a fifth GBM feature. That is the whole exposure.

Two arms, identical folds and targets, three labellings of the KSSL rows:

  A  as reference rows: non-KSSL targets, profile-grouped folds. Does
     mislabelling 12 % of the reference change everyone else's answer?
  B  as targets: KSSL targets with the whole KSSL source held out, so only
     the target's OWN declared type acts. This is the user-facing case --
     someone with a sieved dry end saying "undisturbed" or "disturbed".

No rows are added or dropped in any arm.

Result (2026-10, 480 non-KSSL and 410 KSSL targets, n_mc=30): nothing
significant, and the arms disagree in sign -- relabelling "disturbed" gives
arm A +1.3 pp (p=0.42) and arm B -2.9 pp (p=0.24); "unknown" gives +0.6
(p=0.73) and -2.7 (p=0.27). The label is not inert, though: 162 of 410 arm B
answers change while accuracy does not, i.e. it shuffles predictions without
improving them, consistent with verify_sample_type.py scoring the feature at
+0.0 pp (p=1.00) own laboratory and -0.9 (p=0.58) unseen. DECISION: leave all
2,530 rows labelled "undisturbed" and document the mixed provenance instead.

Usage:  python verify_kssl_label.py [n_per_class] [n_mc]
"""
import sys

import numpy as np
import pandas as pd
import swcc_texture as st
from verify_common import mcnemar
from verify_groups import GROUP

H = np.concatenate([[0.0], np.logspace(-1, np.log10(1500.0), 25)])
LABELS = ("undisturbed", "disturbed", "unknown")
N_PER_CLASS = int(sys.argv[1]) if len(sys.argv) > 1 else 30
N_MC = int(sys.argv[2]) if len(sys.argv) > 2 else 30
N_FOLDS = 5


def curve(r):
    return st.vg_theta(H, r.thetar, r.thetas, r.alpha_kpa, r.n)


df = st.load_reference_df()
df["source_db"] = df.source_db.fillna("KSSL")
df["profile_id"] = df.profile_id.fillna("solo_" + df.layer_id.astype(str))
df["sample_type"] = df.sample_type.fillna("unknown")
is_k = (df.source_db == "KSSL").to_numpy()
print(f"merged {len(df)} layers, KSSL {is_k.sum()}")
nk_ks = df.ksat_cmh[is_k].notna() & (df.ksat_cmh[is_k] > 0)
assert not nk_ks.any(), "KSSL has a measured Ks -- the Ks arm is NOT vacuous"
print("confirmed: no KSSL row carries a measured Ks, so Ks cannot be affected\n")

rng = np.random.default_rng(0)
# Arm A targets: non-KSSL soils.
pool_a = df[~is_k]
tg_a = pd.concat([g.sample(min(len(g), N_PER_CLASS), random_state=0)
                  for _, g in pool_a.groupby("texture_class")]).reset_index(drop=True)
# Arm B targets: KSSL soils.
tg_b = pd.concat([g.sample(min(len(g), N_PER_CLASS), random_state=0)
                  for _, g in df[is_k].groupby("texture_class")]).reset_index(drop=True)
print(f"arm A targets {len(tg_a)} (non-KSSL), arm B targets {len(tg_b)} (KSSL); "
      f"n_mc={N_MC}\n")

profs = np.array(sorted(tg_a.profile_id.unique())); rng.shuffle(profs)
folds_a = [set(x) for x in np.array_split(profs, N_FOLDS)]

res = {}
for lab in LABELS:
    d = df.copy()
    d.loc[is_k, "sample_type"] = lab

    # Arm A: profile-grouped, KSSL always in the reference.
    pred = np.empty(len(tg_a), dtype=object)
    for held in folds_a:
        sel = np.where(tg_a.profile_id.isin(held).to_numpy())[0]
        if not len(sel):
            continue
        train = d[~d.profile_id.isin(held)]
        ref = st.GshpReference(df=train.reset_index(drop=True))
        clf = st.TextureGBM(df=train.reset_index(drop=True))
        for j in sel:
            r = tg_a.iloc[j]
            pred[j] = st.estimate(H, curve(r), ref=ref, n_mc=N_MC, clf=clf,
                                  sample_type=r.sample_type)["texture_class"]
    # Arm B: whole KSSL source held out; only the target's own label acts.
    train_b = d[~is_k]
    ref_b = st.GshpReference(df=train_b.reset_index(drop=True))
    clf_b = st.TextureGBM(df=train_b.reset_index(drop=True))
    pred_b = np.array([st.estimate(H, curve(r), ref=ref_b, n_mc=N_MC, clf=clf_b,
                                   sample_type=lab)["texture_class"]
                       for r in tg_b.itertuples()], dtype=object)
    res[lab] = (pred, pred_b)
    ta, tb = tg_a.texture_class.to_numpy(), tg_b.texture_class.to_numpy()
    ga = np.mean([GROUP[p] == GROUP[t] for p, t in zip(pred, ta)])
    gb = np.mean([GROUP[p] == GROUP[t] for p, t in zip(pred_b, tb)])
    print(f"KSSL labelled {lab:<12s}  arm A exact {np.mean(pred == ta)*100:5.1f} % "
          f"group {ga*100:5.1f} %   |   arm B exact {np.mean(pred_b == tb)*100:5.1f} % "
          f"group {gb*100:5.1f} %", flush=True)

print("\nMcNemar against the status quo (undisturbed):")
ta, tb = tg_a.texture_class.to_numpy(), tg_b.texture_class.to_numpy()
for lab in LABELS[1:]:
    for arm, truth, i in (("A", ta, 0), ("B", tb, 1)):
        base = res["undisturbed"][i] == truth
        alt = res[lab][i] == truth
        print(f"  arm {arm}, {lab:<12s} {(alt.mean()-base.mean())*100:+5.1f} pp "
              f"p={mcnemar(base, alt):.3f}   "
              f"(answers differing: {int((res['undisturbed'][i] != res[lab][i]).sum())}"
              f"/{len(truth)})")
