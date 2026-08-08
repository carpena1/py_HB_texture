"""One-time preprocessing of the UNSODA 2.0 database.

Reads the UNSODA Access database (data/unsoda_raw/unsoda.mdb, Nemes et al.
2001, https://doi.org/10.1016/S0022-1694(01)00465-6) and writes a per-soil
reference table data/unsoda_reference.csv in the SAME schema as
data/gshp_reference.csv, so the two can be concatenated.

UNSODA ships raw retention data rather than fitted parameters, so the van
Genuchten parameters here are fitted with the same routine the tool uses
(Mualem m = 1 - 1/n).

Units in UNSODA (verified empirically in main(), see notes there):
    preshead        cm of water   -> h [kPa]   = preshead * 0.0980665
    k_sat           cm/day        -> ksat_cmh  = k_sat / 24
    particle_size   micrometres, particle_fraction cumulative (0-1)
    depth_upper/lower  cm
Requires mdbtools (`brew install mdbtools`) to export the .mdb tables.
"""

import io
import subprocess

import numpy as np
import pandas as pd

import swcc_texture as st

MDB = "data/unsoda_raw/unsoda.mdb"
CM_TO_KPA = 0.0980665
MIN_POINTS = 5          # need >= 5 points to fit 4 parameters
MAX_RMSE = 0.03         # reject poor fits (m3/m3)


def mdb_table(name):
    out = subprocess.run(["mdb-export", MDB, name], capture_output=True,
                         text=True, check=True).stdout
    return pd.read_csv(io.StringIO(out), low_memory=False)


def particle_fractions(psd):
    """Interpolate each soil's cumulative particle-size distribution at the
    USDA cutoffs: clay < 2 um, silt 2-50 um, sand > 50 um. Returns % values."""
    rows = {}
    for code, g in psd.groupby("code"):
        g = g.dropna(subset=["particle_size", "particle_fraction"])
        g = g.sort_values("particle_size")
        if len(g) < 3:
            continue
        x = g["particle_size"].to_numpy(float)
        y = g["particle_fraction"].to_numpy(float)
        if y.max() > 1.5:            # some entries are in %
            y = y / 100.0
        if x.min() > 2 or x.max() < 50:
            continue                 # cannot bracket both cutoffs
        clay = np.interp(2.0, x, y)
        silt50 = np.interp(50.0, x, y)
        sand = 1.0 - silt50
        silt = silt50 - clay
        if min(clay, silt, sand) < -0.01:
            continue
        tot = clay + silt + sand
        rows[code] = (100 * sand / tot, 100 * silt / tot, 100 * clay / tot)
    return rows


def main():
    gen = mdb_table("general")
    props = mdb_table("soil_properties")
    ht = mdb_table("lab_drying_h-t")
    psd = mdb_table("particle_size")

    fracs = particle_fractions(psd)
    gen = gen.set_index("code")
    props = props.set_index("code")

    records = []
    n_fit_fail = 0
    for code, g in ht.groupby("code"):
        g = g.dropna(subset=["preshead", "theta"])
        g = g[(g.preshead >= 0) & (g.theta > 0)]
        if len(g) < MIN_POINTS or code not in fracs or code not in gen.index:
            continue
        h = g["preshead"].to_numpy(float) * CM_TO_KPA
        theta = g["theta"].to_numpy(float)
        try:
            popt, _ = st.fit_vg(h, theta)
        except Exception:
            n_fit_fail += 1
            continue
        thetar, thetas, la, ln1 = popt
        alpha, n = 10.0 ** la, 1.0 + 10.0 ** ln1
        rmse = float(np.sqrt(np.mean(
            (st.vg_theta(h, thetar, thetas, alpha, n) - theta) ** 2)))
        if not np.isfinite(rmse) or rmse > MAX_RMSE or n <= 1.0:
            n_fit_fail += 1
            continue

        sand, silt, clay = fracs[code]
        # Use the class implied by the measured fractions, so the label is
        # consistent with the USDA definition used everywhere else.
        cls = st.usda_class(sand, silt, clay)

        ksat = props["k_sat"].get(code, np.nan)
        ksat_cmh = ksat / 24.0 if pd.notna(ksat) and ksat > 0 else np.nan

        du = gen["depth_upper"].get(code, np.nan)
        dl = gen["depth_lower"].get(code, np.nan)
        depth = (du + dl) / 2.0 if pd.notna(du) and pd.notna(dl) else np.nan

        # UNSODA soil codes are NNNL, where the leading digits identify the
        # profile/site and the last digit the horizon within it (e.g. 1010-1015
        # are successive horizons of one Troup profile).
        records.append(dict(layer_id=f"UNSODA_{code}",
                            profile_id=f"UNSODA_P{code // 10}",
                            texture_class=cls,
                            alpha_kpa=alpha, n=n, thetar=thetar, thetas=thetas,
                            sand=sand, silt=silt, clay=clay,
                            ksat_cmh=ksat_cmh, depth_cm=depth, rmse=rmse))

    out = pd.DataFrame(records)
    out.to_csv("data/unsoda_reference.csv", index=False)

    print(f"UNSODA soils with usable curves: {len(out)}  (rejected fits: {n_fit_fail})")
    print(f"  with ksat:  {out['ksat_cmh'].notna().sum()}")
    print(f"  with depth: {out['depth_cm'].notna().sum()}")
    print(f"  median fit RMSE: {out['rmse'].median():.4f}")
    print(out["texture_class"].value_counts())

    # --- unit sanity check -------------------------------------------------
    # 1) Do PSD-derived classes agree with UNSODA's own stated texture names?
    stated = gen["texture"].astype(str).str.strip().str.lower()
    codes = [int(s.split("_")[1]) for s in out["layer_id"]]
    agree = [stated.get(c) == cls for c, cls in zip(codes, out["texture_class"])]
    print(f"\nPSD-derived class matches UNSODA's stated texture for "
          f"{sum(agree)}/{len(agree)} soils "
          f"(confirms particle_size is um and fractions are cumulative)")

    # 2) Are fitted alphas on the same scale as GSHP's (confirms cm -> kPa)?
    gshp = pd.read_csv("data/gshp_reference.csv")
    print("\nclass-median alpha [kPa^-1] and n, UNSODA vs GSHP:")
    for cls in ["sand", "sandy loam", "silt loam", "clay"]:
        u = out[out.texture_class == cls]
        gg = gshp[gshp.texture_class == cls]
        if len(u) and len(gg):
            print(f"  {cls:<12s} alpha {u.alpha_kpa.median():7.3f} vs "
                  f"{gg.alpha_kpa.median():7.3f}   n {u.n.median():5.2f} vs "
                  f"{gg.n.median():5.2f}   (n_soils {len(u)} vs {len(gg)})")


if __name__ == "__main__":
    main()
