"""Should andic properties join the neighbour search (fractions, Ks, the kNN
class vote)?

Same targets and folds as verify_volcanic.py: every volcanic layer (known or
probable), a control draw of other soils, and the Canary set as a new
source. Each target is matched on its fitted parameters (no Monte Carlo, so
the arms differ only in how neighbours are chosen):

  base          the shipped search: four vG parameters
  dim 1, dim 2  + a distance penalty of lambda standard units between a
                target and a reference soil that differ in andic
                properties (reference yes/likely = 1, everything else 0;
                the target states its own value)
  filter        andic targets take neighbours only from andic reference
                soils; other targets only from the rest

Fractions and the class vote use every neighbour; Ks uses the undisturbed
ones for an undisturbed target, as the tool does.

Usage:  python verify_andic_knn.py [n_per_class] [seed]
"""

import sys

import numpy as np
import pandas as pd

import ksat_metrics as km
import swcc_texture as st
from verify_common import GROUP, mcnemar
from scipy.stats import wilcoxon

K = 30
N_FOLDS = 5
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 0
ARMS = {"base": ("dim", 0.0), "dim 1": ("dim", 1.0), "dim 2": ("dim", 2.0),
        "filter": ("filter", None)}


def andic_code(df):
    return df.andic.isin(["yes", "likely"]).to_numpy(float)


def knn(ref, a_ref, t, a_t, excluded, arm):
    """Class vote, weighted fraction mean and Ks median for one target."""
    mode, lam = ARMS[arm]
    f = (np.array([np.log10(t.alpha_kpa), np.log10(t.n - 1.0), t.thetar,
                   t.thetas]) - ref.mean) / ref.std
    d2 = ((ref.z - f) ** 2).sum(axis=1)
    if mode == "dim":
        d2 = d2 + (lam * (a_ref - a_t)) ** 2
    d = np.sqrt(d2)
    d = np.where(excluded, np.inf, d)
    if mode == "filter":
        d = np.where(a_ref == a_t, d, np.inf)

    def pick(dist):
        idx = np.argpartition(dist, K)[:K]
        idx = idx[np.isfinite(dist[idx])]
        w = 1.0 / (dist[idx] ** 2 + 1e-6) * ref.class_weight[idx] * ref.row_weight[idx]
        return idx, w / w.sum()

    idx, w = pick(d)
    votes = pd.Series(w).groupby(ref.classes[idx]).sum()
    frac = (ref.fractions[idx] * w[:, None]).sum(axis=0)
    if t.sample_type == "undisturbed":
        idx, w = pick(np.where(ref.sample_type == "undisturbed", d, np.inf))
    ok = np.isfinite(ref.ksat[idx])
    ks = np.nan
    if ok.any():
        lk, ww = np.log10(ref.ksat[idx][ok]), w[ok]
        o = np.argsort(lk)
        ks = 10 ** np.interp(0.5 * ww.sum(), np.cumsum(ww[o]), lk[o])
    return votes.idxmax(), frac, ks


def run(ref_df, ref, a_ref, tg, key, extra):
    n = len(tg) + (len(extra) if extra is not None else 0)
    out = {a: dict(cls=np.empty(n, object), frac=np.full((n, 3), np.nan),
                   ks=np.full(n, np.nan)) for a in ARMS}
    allt = pd.concat([tg, extra]) if extra is not None else tg
    a_all = andic_code(allt)
    groups = np.array(sorted(tg[key].unique()))
    np.random.default_rng(SEED).shuffle(groups)
    chunks = list(np.array_split(groups, N_FOLDS))
    for j in range(n):
        t = allt.iloc[j]
        if j < len(tg):
            chunk = next(c for c in chunks if t[key] in set(c))
            excl = ref_df[key].isin(chunk).to_numpy()
        else:
            excl = np.zeros(len(ref_df), bool)
        for arm in ARMS:
            c, f, k = knn(ref, a_ref, t, a_all[j], excl, arm)
            out[arm]["cls"][j], out[arm]["frac"][j], out[arm]["ks"][j] = c, f, k
    return allt, out


def report(title, allt, out, mask):
    truth = allt.texture_class.to_numpy()[mask]
    meas = allt[["sand", "silt", "clay"]].to_numpy(float)[mask]
    obs = allt.ksat_cmh.to_numpy(float)[mask]
    print(f"\n  {title} (n={mask.sum()}, {np.isfinite(obs).sum()} with Ks)")
    b = out["base"]
    b_ok = b["cls"][mask] == truth
    b_err = np.abs(b["frac"][mask] - meas).mean(axis=1)
    b_ks = np.abs(np.log10(b["ks"][mask]) - np.log10(obs))
    for arm in ARMS:
        o = out[arm]
        ok = o["cls"][mask] == truth
        g = np.mean([GROUP[x] == GROUP[y] for x, y in zip(o["cls"][mask], truth)])
        err = np.abs(o["frac"][mask] - meas).mean(axis=1)
        mae = np.abs(o["frac"][mask] - meas).mean(axis=0)
        ke = np.abs(np.log10(o["ks"][mask]) - np.log10(obs))
        fin = np.isfinite(ke) & np.isfinite(b_ks)
        line = (f"    {arm:<7s} kNN class {ok.mean()*100:5.1f} % / group {g*100:5.1f} %   "
                f"fractions MAE {err.mean():4.1f} (sand {mae[0]:4.1f} silt {mae[1]:4.1f} "
                f"clay {mae[2]:4.1f})")
        if fin.sum() > 5:
            line += f"   Ks med |log err| {np.nanmedian(ke):.2f}"
        if arm != "base":
            line += f"\n            vs base: class {(ok.mean()-b_ok.mean())*100:+.1f} pp p={mcnemar(b_ok, ok):.3f}"
            dm = err - b_err
            if np.any(dm != 0):
                line += f"; fractions {dm.mean():+.2f} pts p={wilcoxon(err, b_err).pvalue:.3f}"
            if fin.sum() > 5 and np.any(ke[fin] != b_ks[fin]):
                line += f"; Ks {np.median(ke[fin]) - np.median(b_ks[fin]):+.2f} dex p={wilcoxon(ke[fin], b_ks[fin]).pvalue:.3f}"
        print(line)


def main():
    n_per_class = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    ref_df = st.load_reference_df()
    ref_df["source_db"] = ref_df.source_db.fillna(ref_df.layer_id.astype(str).str.split("_").str[0])
    ref_df["profile_id"] = ref_df.profile_id.fillna("solo_" + ref_df.layer_id.astype(str))
    ref_df = ref_df.reset_index(drop=True)
    ref = st.GshpReference(df=ref_df)
    a_ref = andic_code(ref_df)
    volc = ref_df.volcanic.isin(["known", "probable"])
    rest = ref_df[~volc & (ref_df.volcanic == "unlikely")]
    ctrl = pd.concat([g.sample(min(len(g), n_per_class), random_state=SEED)
                      for _, g in rest.groupby("texture_class")])
    tg = pd.concat([ref_df[volc], ctrl]).reset_index(drop=True)
    try:
        extra = st.load_reference_df("armas").reset_index(drop=True)
    except FileNotFoundError:
        extra = None
    print(f"reference {len(ref_df)}; andic soils in it {int(a_ref.sum())}; "
          f"targets {volc.sum()} volcanic + {len(ctrl)} controls"
          + (f" + {len(extra)} Canary" if extra is not None else ""))
    for key, name in (("profile_id", "profile folds (source in the reference)"),
                      ("source_db", "source folds (new source)")):
        allt, out = run(ref_df, ref, a_ref, tg, key, extra if key == "source_db" else None)
        v = allt.volcanic.isin(["known", "probable"]).to_numpy()
        a = allt.andic.isin(["yes", "likely"]).to_numpy()
        print(f"\n=== {name} ===")
        report("andic targets", allt, out, a)
        report("volcanic, not andic", allt, out, v & ~a)
        report("control (not volcanic)", allt, out, ~v)
        if key == "source_db" and extra is not None:
            c = np.r_[np.zeros(len(tg), bool), np.ones(len(extra), bool)]
            report("Canary (not in the reference)", allt, out, c)


if __name__ == "__main__":
    main()
