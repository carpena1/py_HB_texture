"""Can soils with only two or three retention points still teach the tool?

A curve fit needs five points, so every source that reports water content at
just 0, 33 and 1500 kPa is thrown away -- Canada's National Pedon Database,
EU-HYDI's Spanish and Greek contributors, and the samples other contributors
report too thinly to fit.

With water contents at fixed heads as predictors (swcc_texture.HEADS_CM,
verify_alt_predictors.py) those soils become usable: three of the eight heads
are saturation, 330 cm (33 kPa) and 15000 cm (1500 kPa), and the
gradient-boosted classifier takes missing features natively. They cannot join
the neighbour search, which needs a complete feature vector, so they are
classifier-only training rows.

This script measures whether adding them helps. Arms, on identical targets
and folds:

  reference        the merged reference alone, heads from each layer's vG fit
  + partial        the same plus partial-head rows, which carry a texture
                   class and whatever heads their points give

Partial rows come from:
  npdb     Canada's NPDB (Open Government Licence), Physical.csv
  euhydi   EU-HYDI samples its preparation drops for too few points, which
           is where Spain and Greece sit (restricted; local only)

A head takes a measured point within TOL log units of it, or log-linear
interpolation between the two points bracketing it; heads outside a sample's
measured range stay missing.

Usage:  python verify_partial_heads.py [profile|source] [n_per_class]
"""

import os
import sys

import numpy as np
import pandas as pd

import swcc_texture as st
from verify_common import GROUP, mcnemar

NPDB = os.path.join("data", "npdb_raw", "Physical.csv")
CM_TO_KPA = 0.0980665
TOL = 0.1                  # log10 units: a point this close to a head is it
N_FOLDS = 5
HEAD_COLS = [f"th{int(c)}" for c in st.HEADS_CM]


def heads_from_points(h_kpa, theta):
    """Water content at HEADS_CM from measured points, NaN where unreachable."""
    out = np.full(len(st.HEADS_CM), np.nan)
    h = np.asarray(h_kpa, float)
    t = np.asarray(theta, float)
    ok = np.isfinite(h) & np.isfinite(t)
    h, t = h[ok], t[ok]
    if not len(h):
        return out
    for j, head in enumerate(st.HEADS_KPA):
        if head == 0.0:
            if (h == 0).any():
                out[j] = t[h == 0].max()
            continue
        pos = h > 0
        if not pos.any():
            continue
        lh, lt = np.log10(h[pos]), t[pos]
        o = np.argsort(lh)
        lh, lt = lh[o], lt[o]
        target = np.log10(head)
        near = np.abs(lh - target).argmin()
        if abs(lh[near] - target) <= TOL:
            out[j] = lt[near]
        elif lh[0] < target < lh[-1]:
            out[j] = np.interp(target, lh, lt)
    return out


def npdb_partial():
    """Canadian NPDB rows with a texture and at least one retention value."""
    if not os.path.exists(NPDB):
        print(f"  {NPDB} not found -- skipping Canada", file=sys.stderr)
        return pd.DataFrame()
    p = pd.read_csv(NPDB, low_memory=False)
    num = lambda c: pd.to_numeric(p[c], errors="coerce")
    sand, silt, clay = num("T_SAND"), num("T_SILT"), num("T_CLAY")
    tot = sand + silt + clay
    keep = tot.between(95, 105) & (sand >= 0) & (clay >= 0)
    # volumetric %, the spec's "percent by total soil volume"
    vals = {0.0: num("RETN_0KP"), 330.0: num("RETN_33KP"),
            15000.0: num("RETN_1500K")}
    rows = []
    for i in np.where(keep)[0]:
        th = {c: v.iloc[i] / 100.0 for c, v in vals.items()
              if np.isfinite(v.iloc[i]) and v.iloc[i] > 0}
        if not th:
            continue
        s, si, c = (100 * x.iloc[i] / tot.iloc[i] for x in (sand, silt, clay))
        r = dict(texture_class=st.usda_class(s, si, c), source_db="NPDB_Canada")
        r.update({f"th{int(k)}": np.nan for k in st.HEADS_CM})
        r.update({f"th{int(k)}": v for k, v in th.items()})
        rows.append(r)
    return pd.DataFrame(rows)


def euhydi_partial():
    """EU-HYDI samples dropped for too few points (Spain, Greece, others)."""
    try:
        import prepare_euhydi as pe
    except ImportError:
        return pd.DataFrame()
    if not os.path.exists(pe.RAW):
        print("  EU-HYDI raw database not here -- skipping", file=sys.stderr)
        return pd.DataFrame()
    have = set()
    if os.path.exists(pe.OUT):
        have = set(pd.read_csv(pe.OUT).layer_id.astype(str))
    ret = pe.table("RET")
    for c in ("HEAD", "THETA"):
        ret[c] = pd.to_numeric(ret[c], errors="coerce")
    ret = ret[ret.FLAG.astype(str) != "False"]
    ret = ret[ret.HEAD.notna() & ret.THETA.notna() & (ret.HEAD >= 0)]
    psd = pe.table("PSD_EST").set_index("SAMPLE_ID")
    psize = pe.table("PSIZE")
    for c in ("P_SIZE", "P_PERCENT"):
        psize[c] = pd.to_numeric(psize[c], errors="coerce")
    classes = pe.usda_from_size_classes(psize)
    gen = pe.table("GENERAL").set_index("PROFILE_ID")
    rows = []
    for sid, g in ret.groupby("SAMPLE_ID"):
        if f"EUHYDI_{sid}" in have:
            continue
        frac = None
        if sid in psd.index:
            q = psd.loc[sid]
            q = q.iloc[0] if isinstance(q, pd.DataFrame) else q
            if {q.USCLAY_CODE, q.USSILT_CODE, q.USSAND_CODE} <= pe.PSD_OK:
                frac = (float(q.USSAND), float(q.USSILT), float(q.USCLAY))
        if (frac is None or not np.isfinite(sum(frac))) and sid in classes:
            frac = classes[sid]
        if frac is None or not np.isfinite(sum(frac)) \
                or not (95 <= sum(frac) <= 105):
            continue
        s, si, c = (100 * x / sum(frac) for x in frac)
        th = heads_from_points(g.HEAD.to_numpy(float) * CM_TO_KPA,
                               g.THETA.to_numpy(float))
        if not np.isfinite(th).any():
            continue
        pid = g.PROFILE_ID.iloc[0]
        gi = gen.loc[pid] if pid in gen.index else None
        src = (gi.SOURCE if gi is not None and pd.notna(gi.SOURCE)
               else g.SOURCE.iloc[0])
        r = dict(texture_class=st.usda_class(s, si, c),
                 source_db=f"EUHYDI_{src}")
        r.update(dict(zip(HEAD_COLS, th)))
        rows.append(r)
    return pd.DataFrame(rows)


def head_frame(df):
    """Reference layers as head features, from each layer's vG fit."""
    x = st._head_features(df.thetar.to_numpy(), df.thetas.to_numpy(),
                          df.alpha_kpa.to_numpy(), df.n.to_numpy(),
                          st.HEADS_KPA)
    out = pd.DataFrame(x, columns=HEAD_COLS, index=df.index)
    out["texture_class"] = df.texture_class.to_numpy()
    out["source_db"] = df.source_db.to_numpy()
    return out


def gbm(x, y):
    from sklearn.ensemble import HistGradientBoostingClassifier
    return HistGradientBoostingClassifier(
        max_iter=400, learning_rate=0.06, max_leaf_nodes=31,
        l2_regularization=1.0, early_stopping=True, validation_fraction=0.15,
        random_state=0, class_weight="balanced").fit(x, y)


def main():
    design = sys.argv[1] if len(sys.argv) > 1 else "profile"
    n_per_class = int(sys.argv[2]) if len(sys.argv) > 2 else 150

    g = st.load_reference_df("merged").reset_index(drop=True)
    g["source_db"] = g.source_db.fillna("KSSL")
    g["profile_id"] = g.profile_id.fillna("solo_" + g.layer_id.astype(str))
    parts = []
    for name, f in (("Canada NPDB", npdb_partial), ("EU-HYDI dropped", euhydi_partial)):
        p = f()
        if len(p):
            n_heads = np.isfinite(p[HEAD_COLS].to_numpy(float)).sum(axis=1)
            print(f"{name}: {len(p)} partial rows, "
                  f"{np.mean(n_heads):.1f} heads each "
                  f"(>=3 heads: {(n_heads >= 3).sum()})")
            parts.append(p)
    partial = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    if partial.empty:
        raise SystemExit("no partial-head rows available")
    print(f"partial rows in total: {len(partial)}")
    print(partial.texture_class.value_counts().to_string(), "\n")

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

    ref_x = head_frame(g)
    tgt_x = head_frame(tg)
    out = {a: np.empty(len(tg), object) for a in ("reference", "+ partial")}
    for key, held in fold_list:
        sel = np.where(tg[key].isin(held))[0]
        if not len(sel):
            continue
        train = ref_x[~g[key].isin(held).to_numpy()]
        # a held-out laboratory is held out of the partial rows too
        pt = partial[~partial.source_db.isin(held)] if key == "source_db" \
            else partial
        for arm, tr in (("reference", train),
                        ("+ partial", pd.concat([train, pt], ignore_index=True))):
            m = gbm(tr[HEAD_COLS].to_numpy(float),
                    tr.texture_class.to_numpy())
            out[arm][sel] = m.predict(tgt_x.iloc[sel][HEAD_COLS].to_numpy(float))
        print(f"  fold {sorted(held)[0]}: {len(sel)} targets", flush=True)

    truth = tg.texture_class.to_numpy()
    print(f"\n=== {design} design, {len(tg)} targets ===")
    base = None
    for arm in ("reference", "+ partial"):
        ex = out[arm] == truth
        gp = np.array([GROUP[a] == GROUP[b] for a, b in zip(out[arm], truth)])
        line = f"{arm:12s} exact {ex.mean() * 100:5.1f}  group {gp.mean() * 100:5.1f}"
        if base is not None:
            line += (f"   {(ex.mean() - base[0].mean()) * 100:+5.1f} pp "
                     f"(p={mcnemar(base[0], ex):.3f})  group "
                     f"{(gp.mean() - base[1].mean()) * 100:+5.1f} "
                     f"(p={mcnemar(base[1], gp):.3f})")
        print(line)
        base = (ex, gp) if arm == "reference" else base
    print("\nper class (exact recall, reference -> + partial):")
    for c in st.USDA_CLASSES:
        m = truth == c
        if m.sum() < 10:
            continue
        a = (out["reference"][m] == c).mean() * 100
        b = (out["+ partial"][m] == c).mean() * 100
        extra = int((partial.texture_class == c).sum())
        print(f"  {c:18s} n={m.sum():4d}  {a:5.1f} -> {b:5.1f} "
              f"({b - a:+5.1f})   partial rows added: {extra}")


if __name__ == "__main__":
    main()
