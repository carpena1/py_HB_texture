"""Does region-matched reference data actually improve classification?

This estimates what HYPRES would buy us *before* obtaining it. HYPRES's value
proposition is that it would add ~4,486 European horizons to a reference in
which Europe is thinly represented (778 of 9,996 curated GSHP layers, 8 %).

We cannot test HYPRES directly without the data, but we can test the mechanism
by removing European data from the reference and measuring the damage to
European target soils. Three conditions, all with the target's own profile
excluded, compared pairwise on the same targets:

    A  full        reference as-is                  (region match present)
    B  no-Europe   every European layer removed     (region match lost,
                                                     reference 778 smaller)
    C  size-match  778 random NON-European layers   (same size loss as B,
                   removed                           region match kept)

B vs C is the informative contrast: it separates "losing European soils" from
"having a smaller reference". If B is clearly worse than C, region-matched
data carries real signal and adding HYPRES's European horizons should help
European soils. If B ~ C, HYPRES would mostly add bulk, not value.

A non-European control group is also run: removing European data should NOT
hurt non-European soils. If it does, the effect is not regional.

Usage:  python verify_region.py [n_mc]
"""

import sys

import numpy as np
import pandas as pd

import swcc_texture as st
from verify_groups import GROUP

H = np.arange(0.0, 150.0 + 0.1, 5.0)
EU_BBOX = dict(lat=(34, 72), lon=(-25, 45))


def europe_mask(df):
    return (df.lat.between(*EU_BBOX["lat"]) & df.lon.between(*EU_BBOX["lon"])
            ).fillna(False)


def predict(ref_df, targets, n_mc):
    """Predict every target with its own profile excluded from ref_df."""
    ref = st.GshpReference(df=ref_df.reset_index(drop=True))
    preds = []
    for row in targets.itertuples():
        theta = st.vg_theta(H, row.thetar, row.thetas, row.alpha_kpa, row.n)
        ref.set_excluded_profile(row.profile_id)
        preds.append(st.estimate(H, theta, ref=ref,
                                 n_mc=n_mc)["texture_class"])
    return np.array(preds)


def mcnemar(a_ok, b_ok):
    from math import comb
    b = int(np.sum(a_ok & ~b_ok))
    c = int(np.sum(~a_ok & b_ok))
    n = b + c
    if n == 0:
        return b, c, 1.0
    lo = min(b, c)
    return b, c, min(1.0, sum(comb(n, i) for i in range(lo + 1)) / 2 ** n * 2)


def evaluate(preds, truth):
    cls = preds == truth
    grp = np.array([GROUP[p] == GROUP[t] for p, t in zip(preds, truth)])
    return cls, grp


def main():
    n_mc = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    rng = np.random.default_rng(0)

    df = pd.read_csv(st.REFERENCE_CSV)
    eu = europe_mask(df)
    print(f"curated reference: {len(df)} layers, {eu.sum()} European "
          f"({100*eu.mean():.1f} %), {df[eu].profile_id.nunique()} EU profiles")

    eu_targets = df[eu]
    # Non-European control group, matched in size to the EU target set.
    non_eu_pool = df[~eu]
    ctrl_targets = non_eu_pool.sample(len(eu_targets), random_state=1)

    # Condition C: drop as many random non-European layers as B drops.
    drop_n = int(eu.sum())
    drop_idx = rng.choice(non_eu_pool.index.to_numpy(), size=drop_n,
                          replace=False)

    conditions = {
        "A full": df,
        "B no-Europe": df[~eu],
        "C size-matched": df.drop(index=drop_idx),
    }

    for group_name, targets in [("European targets", eu_targets),
                                ("non-European control", ctrl_targets)]:
        print(f"\n{'=' * 64}\n{group_name}  (n = {len(targets)})\n{'=' * 64}")
        truth = targets.texture_class.to_numpy()
        res = {}
        for cname, cdf in conditions.items():
            preds = predict(cdf, targets, n_mc)
            cls, grp = evaluate(preds, truth)
            res[cname] = (cls, grp)
            print(f"  {cname:<16s} ref={len(cdf):5d}  "
                  f"class {cls.mean()*100:5.1f} %   group {grp.mean()*100:5.1f} %")

        print("\n  paired contrasts:")
        pairs = [("A full", "B no-Europe", "cost of losing Europe (+size)"),
                 ("C size-matched", "B no-Europe", "cost of losing Europe "
                                                   "(size held equal)"),
                 ("A full", "C size-matched", "cost of a smaller reference")]
        for x, y, label in pairs:
            for lvl, i in [("class", 0), ("group", 1)]:
                xa, yb = res[x][i], res[y][i]
                b, c, p = mcnemar(xa, yb)
                print(f"    {label:<34s} {lvl:<6s} "
                      f"{(yb.mean()-xa.mean())*100:+6.1f}pp  p={p:6.3f}")


if __name__ == "__main__":
    main()
