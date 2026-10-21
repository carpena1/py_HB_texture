"""Do the former-Soviet-school particle-size data read as too fine?

A. Nemes (email, 2026-10-09): laboratories working to the former Russian
school's standards sieved at 1 mm, not 2 mm, so their "sand" is the
50-1000 um fraction. Reported on a <1 mm basis and renormalised to 100 %,
that makes sand too low and silt and clay too high by the same token. He
notes the data were never converted, in EU-HYDI or in the global
compilations, and suggests looking for the signal directly: do those soils'
curves look like those of sandier soils elsewhere?

The tool answers exactly that question, because it reads fractions from the
curve. If the reported sand of a source is too low, its curves should look
sandier than its own labels say, i.e.

    residual = predicted sand - measured sand

should be more positive for that source than for comparable laboratories.

Three groups:
  suspect    Poland (Lamorski) and Russia (Shein, and Russia_EGRPR, which
             reaches the reference through GSHP and is not in EU-HYDI)
  ambiguous  Hungary, Czechia, Slovakia -- the same school historically,
             but not named in the email
  control    the other European laboratories in the reference

Every target is predicted with its WHOLE SOURCE hidden, so a biased source
cannot be matched against its own rows (the project's rule: profile-level
LOO cannot clear a source-level effect). The residual depends on where a
soil sits in texture space -- sandy soils have their sand under-predicted
and clayey soils over-predicted -- so the control group is also reweighted
to the suspect group's measured-sand distribution before the groups are
compared. Confidence intervals bootstrap over SOURCES, not rows.

Usage:  python verify_psd_standard.py [n_per_source] [n_mc]
"""

import sys

import numpy as np
import pandas as pd

import swcc_texture as st
from verify_common import COLS

MIN_SOURCE = 80
H = np.arange(0.0, 150.0 + 0.1, 5.0)
REF_COLS = COLS + ["source_db"]
SUSPECT = ["EUHYDI_Lamorski", "EUHYDI_Shein", "Russia_EGRPR"]
AMBIGUOUS = ["EUHYDI_Mako", "EUHYDI_Matula", "EUHYDI_Houskova"]
# European laboratories measuring to the 2 mm (USDA/FAO) standard.
CONTROL_EXTRA = ["Hohenbrink", "ZALF_database", "Belgium_database",
                 "Swiss_database"]
SAND_BINS = np.arange(0.0, 100.0 + 1e-9, 10.0)
N_BOOT = 2000


def predict_sand_clay(ref, targets, mask, n_mc):
    """Fractions for each target with its whole source hidden."""
    sand, clay = [], []
    for row in targets.itertuples():
        theta = st.vg_theta(H, row.thetar, row.thetas, row.alpha_kpa, row.n)
        ref.set_excluded(mask)
        f = st.estimate(H, theta, ref=ref, n_mc=n_mc)["fractions"]
        sand.append(f["sand"]); clay.append(f["clay"])
    return np.array(sand, float), np.array(clay, float)


def matched_mean(res, sand, weights_from):
    """Mean residual after reweighting to another group's sand distribution."""
    tgt = np.histogram(weights_from, bins=SAND_BINS)[0].astype(float)
    tgt /= tgt.sum()
    idx = np.clip(np.digitize(sand, SAND_BINS) - 1, 0, len(SAND_BINS) - 2)
    out, wsum = 0.0, 0.0
    for b in range(len(SAND_BINS) - 1):
        m = idx == b
        if m.any() and tgt[b] > 0:
            out += tgt[b] * res[m].mean(); wsum += tgt[b]
    return out / wsum if wsum else np.nan


def boot_ci(by_source, rng):
    """Bootstrap the mean over sources."""
    keys = list(by_source)
    draws = [np.mean([by_source[keys[i]] for i in rng.integers(0, len(keys), len(keys))])
             for _ in range(N_BOOT)]
    return np.percentile(draws, [2.5, 97.5])


def main():
    n_per_source = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    n_mc = int(sys.argv[2]) if len(sys.argv) > 2 else 40

    df = st.load_reference_df()
    df["source_db"] = df.source_db.fillna("KSSL")
    df["profile_id"] = df.profile_id.fillna("solo_" + df.layer_id.astype(str))
    df = df[df.sand.notna() & df.clay.notna()].reset_index(drop=True)
    ref = st.GshpReference(df=df[REF_COLS].reset_index(drop=True))
    src_all = df.source_db.to_numpy()

    sizes = df.source_db.value_counts()
    control = [s for s in sizes.index
               if sizes[s] >= MIN_SOURCE
               and (s.startswith("EUHYDI_") or s in CONTROL_EXTRA)
               and s not in SUSPECT and s not in AMBIGUOUS]
    groups = {"suspect": SUSPECT, "ambiguous": AMBIGUOUS, "control": control}

    print(f"reference {len(df)} layers with fractions; n_mc={n_mc}, "
          f"up to {n_per_source} targets per source\n")
    rows, res = [], {g: {} for g in groups}
    keep = {g: {"res": [], "sand": [], "clay": []} for g in groups}
    for g, sources in groups.items():
        for s in sources:
            if sizes.get(s, 0) < MIN_SOURCE:
                continue
            sub = df[df.source_db == s]
            tg = pd.concat([x.sample(min(len(x), max(1, n_per_source // 12)),
                                     random_state=0)
                            for _, x in sub.groupby("texture_class")])
            ps, pc = predict_sand_clay(ref, tg, src_all == s, n_mc)
            r_sand = ps - tg.sand.to_numpy(float)
            r_clay = pc - tg.clay.to_numpy(float)
            res[g][s] = float(r_sand.mean())
            keep[g]["res"].append(r_sand)
            keep[g]["sand"].append(tg.sand.to_numpy(float))
            keep[g]["clay"].append(r_clay)
            rows.append(dict(group=g, source=s, size=int(sizes[s]), n=len(tg),
                             meas_sand=tg.sand.mean(), pred_sand=ps.mean(),
                             d_sand=r_sand.mean(), d_clay=r_clay.mean()))
            print(f"  {g:<10s} {s:<24s} n={len(tg):4d}  measured sand "
                  f"{tg.sand.mean():5.1f} %  predicted {ps.mean():5.1f} %  "
                  f"residual {r_sand.mean():+5.1f}")

    for g in keep:
        for k in keep[g]:
            keep[g][k] = np.concatenate(keep[g][k]) if keep[g][k] else np.array([])

    print("\n" + "=" * 72)
    print("PER SOURCE (sand, percentage points; + = curve looks sandier than "
          "the label)")
    t = pd.DataFrame(rows)
    print(t.to_string(index=False, float_format=lambda v: f"{v:6.1f}"))

    rng = np.random.default_rng(0)
    print("\nGROUP MEANS (bootstrap over sources, 95 % CI)")
    print(f"{'group':<11s} {'sources':>7s} {'n':>6s} {'sand resid':>11s} "
          f"{'95 % CI':>16s} {'clay resid':>11s}")
    for g in ("suspect", "ambiguous", "control"):
        if not res[g]:
            continue
        lo, hi = boot_ci(res[g], rng)
        # mean over sources, to match the bootstrap; large labs do not count
        # more than small ones.
        print(f"{g:<11s} {len(res[g]):>7d} {len(keep[g]['res']):>6d} "
              f"{np.mean(list(res[g].values())):>+11.1f} "
              f"{f'[{lo:+.1f}, {hi:+.1f}]':>16s} "
              f"{keep[g]['clay'].mean():>+11.1f}")

    print("\nMATCHED ON MEASURED SAND (control reweighted to the suspect "
          "group's sand distribution)")
    for g in ("suspect", "ambiguous"):
        if not res[g]:
            continue
        m_c = matched_mean(keep["control"]["res"], keep["control"]["sand"],
                           keep[g]["sand"])
        m_g = keep[g]["res"].mean()
        print(f"  {g:<10s} {m_g:+.1f} vs control {m_c:+.1f} "
              f"-> difference {m_g - m_c:+.1f} points of sand")

    # A source-level permutation test: is the suspect group's mean residual
    # extreme among all the European sources scored here?
    allsrc = {**res["suspect"], **res["ambiguous"], **res["control"]}
    obs = np.mean([allsrc[s] for s in res["suspect"]])
    vals = np.array(list(allsrc.values()))
    k = len(res["suspect"])
    draws = np.array([rng.choice(vals, k, replace=False).mean()
                      for _ in range(10000)])
    print(f"\npermutation over the {len(vals)} sources: suspect mean "
          f"{obs:+.1f}, p={float((draws >= obs).mean()):.3f} (one-sided)")


if __name__ == "__main__":
    main()
