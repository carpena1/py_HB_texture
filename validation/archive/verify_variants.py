"""Compare reference-database / covariate variants of the tool.

Variants tested (the fitting stage is identical in all of them; only the
reference used for kNN inference changes):

    GSHP                  baseline (data/gshp_reference.csv)
    GSHP+depth            adds sample mid-depth as a 5th kNN feature
    GSHP+UNSODA           merged reference (UNSODA 2.0, vG fitted here)
    GSHP+UNSODA+depth     both

Benchmarks:
    Carsel & Parrish (1988) and ROSETTA class means -- synthetic class-typical
      curves that carry no depth, so only the database variants are scored.
    GSHP leave-one-out and UNSODA leave-one-out -- 12 representative real soils
      each (one per USDA class), with the target soil removed from the
      reference before predicting; these carry a real depth, so all four
      variants are scored.

Reports class-level accuracy and accuracy over the 4 aggregated texture
groups (see verify_groups.py).
"""

import numpy as np
import pandas as pd

import swcc_texture as st
import verify_carsel_parrish as vcp
import verify_rosetta as vr
from verify_groups import GROUP

H = np.arange(0.0, 150.0 + 0.1, 5.0)
COLS = ["layer_id", "profile_id", "texture_class", "alpha_kpa", "n", "thetar",
        "thetas", "sand", "silt", "clay", "ksat_cmh", "depth_cm"]
ORDER = ["sand", "loamy sand", "sandy loam", "loam", "silt", "silt loam",
         "sandy clay loam", "clay loam", "silty clay loam", "sandy clay",
         "silty clay", "clay"]


def feat(d):
    return np.column_stack([np.log10(d.alpha_kpa), np.log10(d.n - 1.0),
                            d.thetar, d.thetas])


def representatives(df, need_depth=True):
    """One soil per class: closest to the class median in vG feature space,
    preferring soils with a measured Ksat (and a depth when needed)."""
    out = {}
    for cls in ORDER:
        g = df[df.texture_class == cls]
        if len(g) == 0:
            continue
        pool = g[g.ksat_cmh.notna()]
        if need_depth:
            p2 = pool[pool.depth_cm.notna()]
            pool = p2 if len(p2) else g[g.depth_cm.notna()]
        if len(pool) == 0:
            pool = g
        med = np.median(feat(g), axis=0)
        z = (feat(pool) - med) / (feat(g).std(axis=0) + 1e-9)
        out[cls] = pool.index[np.argmin((z ** 2).sum(1))]
    return out


def score(preds):
    n = len(preds)
    cls_ok = sum(t == p for t, p in preds)
    grp_ok = sum(GROUP[t] == GROUP[p] for t, p in preds)
    return cls_ok, grp_ok, n


# --- benchmarks -------------------------------------------------------------

def run_carsel(ref_df, use_depth):
    if use_depth:
        return None
    names, curves, _ = vcp.load_benchmark()
    ref = st.GshpReference(df=ref_df)
    return [(n, st.estimate(curves[n]["h"], curves[n]["theta"],
                            ref=ref)["texture_class"]) for n in names]


def run_rosetta(ref_df, use_depth):
    if use_depth:
        return None
    ref = st.GshpReference(df=ref_df)
    out = []
    for n in vr.ORDER:
        tr, ts, a, nn, _ = vr.rosetta_params(n)
        theta = st.vg_theta(H, tr, ts, a, nn)
        out.append((n, st.estimate(H, theta, ref=ref)["texture_class"]))
    return out


def run_loo(ref_df, use_depth, target_df, tag):
    """Leave-one-out over 12 representative soils drawn from `target_df`
    (which is part of `ref_df`); the target is dropped before predicting."""
    reps = representatives(target_df, need_depth=use_depth)
    out = []
    for cls, j in reps.items():
        row = target_df.loc[j]
        theta = st.vg_theta(H, row.thetar, row.thetas, row.alpha_kpa, row.n)
        loo_df = ref_df[ref_df.layer_id != row.layer_id]
        ref = st.GshpReference(df=loo_df.reset_index(drop=True),
                               use_depth=use_depth)
        depth = row.depth_cm if use_depth else None
        out.append((cls, st.estimate(H, theta, ref=ref,
                                     depth=depth)["texture_class"]))
    return out


def main():
    gshp = pd.read_csv(st.REFERENCE_CSV)[COLS]
    unsoda = pd.read_csv("data/unsoda_reference.csv")[COLS]
    merged = pd.concat([gshp, unsoda], ignore_index=True)

    print(f"reference sizes: GSHP {len(gshp)}, UNSODA {len(unsoda)}, "
          f"merged {len(merged)}")
    print(f"  with depth:    GSHP {gshp.depth_cm.notna().sum()}, "
          f"UNSODA {unsoda.depth_cm.notna().sum()}, "
          f"merged {merged.depth_cm.notna().sum()}\n")

    variants = [("GSHP", gshp, False), ("GSHP+depth", gshp, True),
                ("GSHP+UNSODA", merged, False),
                ("GSHP+UNSODA+depth", merged, True)]
    benchmarks = [
        ("Carsel & Parrish", lambda d, u: run_carsel(d, u)),
        ("ROSETTA", lambda d, u: run_rosetta(d, u)),
        ("GSHP LOO", lambda d, u: run_loo(d, u, gshp, "gshp")),
        ("UNSODA LOO", lambda d, u: run_loo(d, u, unsoda, "unsoda")),
    ]

    rows = []
    for vname, vdf, vdepth in variants:
        for bname, fn in benchmarks:
            preds = fn(vdf, vdepth)
            if preds is None:
                rows.append((vname, bname, None, None, None))
                continue
            c, g, n = score(preds)
            rows.append((vname, bname, c, g, n))

    print(f"{'variant':<20s} {'benchmark':<18s} {'class':>10s} {'group':>10s}")
    print("-" * 62)
    for v, b, c, g, n in rows:
        if c is None:
            print(f"{v:<20s} {b:<18s} {'n/a':>10s} {'n/a':>10s}")
        else:
            print(f"{v:<20s} {b:<18s} {f'{c}/{n}':>10s} {f'{g}/{n}':>10s}")

    print("\n(n/a = benchmark carries no depth, so covariate variants are "
          "not scored on it)")

    # Like-for-like: only the two leave-one-out benchmarks can be scored for
    # every variant, so totals are taken over those alone.
    common = {"GSHP LOO", "UNSODA LOO"}
    print("\n=== totals over the two LOO benchmarks (scored for all "
          "variants) ===")
    for vname, _, _ in variants:
        sel = [(c, g, n) for v, b, c, g, n in rows
               if v == vname and b in common and c is not None]
        C, G, N = (sum(x[0] for x in sel), sum(x[1] for x in sel),
                   sum(x[2] for x in sel))
        print(f"  {vname:<20s} class {C}/{N} ({100*C/N:.0f} %)   "
              f"group {G}/{N} ({100*G/N:.0f} %)")


if __name__ == "__main__":
    main()
