"""Build a reference table from the NCSS/KSSL Soil Characterization Database.

Source: USDA-NRCS Kellogg Soil Survey Laboratory, "Lab Data Mart" SQLite
snapshot (https://ncsslabdatamart.sc.egov.usda.gov/database_download.aspx).
US Government work, public domain; the site asks only that the database be
cited. The 5.4 GB snapshot is gitignored -- rebuild data/kssl_reference.csv by
running this script.

WHAT KSSL DOES AND DOES NOT CONTAIN
-----------------------------------
It contains water retention, particle-size analysis and bulk density. It does
NOT contain saturated hydraulic conductivity. lab_analyte does define
"Saturated Conductivity, Replicate 0/1/2" (keys 1803-1805) and "Hydraulic
Conductivity" (key 1913), but their `column_name` is empty and no table in the
snapshot holds the values -- the analytes are declared in the dictionary and
nothing is distributed. So every row written here has ksat_cmh = NaN, and
merging KSSL DILUTES the fraction of the reference that can support a Ksat
estimate. That trade-off is measured, not assumed: see verify_kssl.py.

UNITS
-----
KSSL reports water retention gravimetrically, as percent of oven-dry mass
(a handful of organic layers exceed 700 %, which is impossible volumetrically
and confirms the basis). Volumetric water content is therefore

    theta = (w_gravimetric / 100) * bulk_density

using the 1/3-bar bulk density where available and oven-dry otherwise. Two
independent unit checks run in main(): the PSDA-derived USDA class is compared
against KSSL's own texture_lab field, and the fitted class-median alpha is
compared against GSHP's.

Tensions used (bar -> kPa): 0.06->6, 0.1->10, 1/3->33, 1->100, 2->200,
5->500, 15->1500. The 0-bar and 3-bar columns are skipped: they are populated
for 89 and 0 layers respectively.

alpha is fitted directly against h in kPa here, so it comes out in kPa^-1
natively, matching GSHP (whose published alpha is in m^-1 and is divided by
9.80665 in prepare_gshp.py). Verified: refitting GSHP layers from their raw
(h, theta) points reproduces the published alpha with a median factor of 1.00.

THE SATURATION ANCHOR
---------------------
KSSL's wettest tension is 6 kPa and most layers start at 33 kPa, so nothing
constrains the curve near saturation. Fitting as-is drives thetas up (19 % of
layers came out above 0.75, which is not physical for a mineral soil) and
alpha down by 2-6x relative to GSHP -- a measurement-window artifact, not a
soil difference. Confirmed by truncating GSHP's own measured curves at 6 kPa
and refitting: alpha falls by the same factor and in the same direction
(clay 0.17x, clay loam 0.17x, silty clay loam 0.20x).

The fix is one synthetic point at h = 0 with theta = 0.95 * total porosity,
porosity being 1 - bulk_density_oven_dry / PARTICLE_DENSITY. Validated on
GSHP, where the full-curve answer is known: with the anchor, per-layer
agreement with the true alpha improves from 20 % to 48 % within a factor of
two and the median bias falls from +0.57 to +0.20 dex. It removes most of the
systematic bias but NOT the per-layer scatter, so KSSL rows remain noisier
than GSHP rows.

Two soil-physics caveats behind those constants, both measured on GSHP:

  * theta_s is not porosity. Entrapped air keeps a saturated soil below its
    pore volume; the measured median ratio is 0.946, hence
    SATURATION_FRACTION. (Even so, theta_s exceeds 1 - BD/2.65 in 26 % of
    GSHP layers, which is impossible at that particle density -- see below.)
  * 2.65 g/cm3 is the density of silica and an upper bound for real soil.
    Where GSHP reports porosity and bulk density directly, the implied
    particle density BD/(1-porosity) has median 2.589, and it falls with
    organic matter (corr with organic carbon -0.46; median 1.98 above 5 % OC).
    Organic layers are excluded here by MAX_THETA rather than by modelling
    their particle density, so PARTICLE_DENSITY is kept at the mineral value.
"""

import os
import sqlite3

import numpy as np
import pandas as pd

import swcc_texture as st

DB = "data/kssl_raw/NCSSLabDataMartSQLite.sqlite3"
OUT = "data/kssl_reference.csv"

# column -> matric potential in kPa
RETENTION = {
    "water_retention_6_hundredths": 6.0,
    "water_retention_10th_bar": 10.0,
    "water_retention_third_bar": 33.0,
    "water_retention_1_bar": 100.0,
    "water_retention_2_bar": 200.0,
    "water_retention_5_bar_sieve": 500.0,
    "water_retention_15_bar": 1500.0,
}
MIN_POINTS = 5
MAX_RMSE = 0.03
MAX_THETA = 0.95        # drops organic layers, whose gravimetric water is huge
PARTICLE_DENSITY = 2.65
POROSITY_RANGE = (0.15, 0.90)

# theta_s is NOT equal to total porosity: entrapped air keeps a saturated soil
# below its pore volume. Measured on GSHP (9,958 layers with bulk density),
# theta_s / porosity has median 0.946, so the saturation anchor is placed at
# 0.95 * porosity rather than at porosity itself. Validated by removing the
# wet end from GSHP's own curves and refitting: the median bias in log10 alpha
# falls from +0.278 (anchor = porosity) to +0.204, and agreement within a
# factor of two rises from 45 % to 48 %. Lower fractions reduce the bias
# further (0.90 -> +0.137) but have no physical justification, so 0.95 is used.
SATURATION_FRACTION = 0.95

# KSSL texture_lab codes -> USDA class. Sand-grade qualifiers (very fine,
# fine, medium, coarse, very coarse) do not change the USDA class, so they are
# stripped before lookup. Exact codes win first: 'cl' is clay loam, not
# coarse loam, and 'cosl' is coarse sandy loam, not coarse silt.
TEXTURE_CODES = {
    "s": "sand", "ls": "loamy sand", "sl": "sandy loam", "l": "loam",
    "sil": "silt loam", "si": "silt", "scl": "sandy clay loam",
    "cl": "clay loam", "sicl": "silty clay loam", "sc": "sandy clay",
    "sic": "silty clay", "c": "clay",
}
# Shortest first, so 'cosl' strips 'co' -> 'sl' rather than 'cos' -> 'l'.
SAND_GRADES = ("f", "m", "vf", "vc", "co", "cos")


def lab_code_to_class(code):
    """Map a KSSL texture_lab code such as 'fsl' or 'COSL' to a USDA class."""
    if not isinstance(code, str):
        return None
    c = code.strip().lower()
    if c in TEXTURE_CODES:                       # exact code wins
        return TEXTURE_CODES[c]
    for p in sorted(SAND_GRADES, key=len):       # e.g. fsl -> sl
        if c.startswith(p) and c[len(p):] in TEXTURE_CODES:
            return TEXTURE_CODES[c[len(p):]]
    if c.startswith("l"):                        # lfs, lvfs, lcos -> loamy sand
        for p in sorted(SAND_GRADES, key=len):
            if c[1:] in ("s",) or c[1:] == p + "s":
                return "loamy sand"
    return None


def load():
    q = f"""
    SELECT p.layer_key, p.labsampnum,
           {', '.join('p.' + c for c in RETENTION)},
           p.sand_total, p.silt_total, p.clay_total, p.texture_lab,
           p.bulk_density_third_bar, p.bulk_density_oven_dry,
           l.pedon_key, l.hzn_top, l.hzn_bot, l.hzn_desgn,
           s.latitude_std_decimal_degrees  AS lat,
           s.longitude_std_decimal_degrees AS lon
    FROM lab_physical_properties_vw p
    JOIN lab_layer l ON l.layer_key = p.layer_key
    LEFT JOIN lab_pedon d ON d.pedon_key = l.pedon_key
    LEFT JOIN lab_site  s ON s.site_key  = d.site_key
    """
    with sqlite3.connect(DB) as con:
        df = pd.read_sql(q, con)

    # One layer can appear more than once (different prep codes / sources);
    # keep whichever row carries the most retention points.
    df["npts"] = df[list(RETENTION)].notna().sum(axis=1)
    df = (df.sort_values("npts", ascending=False)
            .drop_duplicates("layer_key", keep="first"))
    return df


def main():
    if not os.path.exists(DB):
        raise SystemExit(f"{DB} not found -- download the Lab Data Mart "
                         f"SQLite snapshot first (see module docstring).")
    df = load()
    print(f"physical-property layers: {len(df)}")

    bd = df.bulk_density_third_bar.fillna(df.bulk_density_oven_dry)
    n_bd33 = int(df.bulk_density_third_bar.notna().sum())
    # Total porosity uses oven-dry bulk density; it supplies the saturation
    # anchor that KSSL's tension range cannot constrain (see module docstring).
    bd_od = df.bulk_density_oven_dry.fillna(df.bulk_density_third_bar)
    porosity = 1.0 - bd_od.to_numpy() / PARTICLE_DENSITY

    rows = []
    skip = dict(points=0, texture=0, bd=0, theta=0, fit=0, rmse=0, porosity=0)
    for row, bdi, por in zip(df.itertuples(), bd.to_numpy(), porosity):
        h, th = [], []
        for col, kpa in RETENTION.items():
            w = getattr(row, col)
            if w is not None and np.isfinite(w):
                h.append(kpa)
                th.append(w)
        if len(h) < MIN_POINTS:
            skip["points"] += 1
            continue
        if not np.isfinite(bdi) or bdi <= 0:
            skip["bd"] += 1
            continue

        h = np.array(h, float)
        theta = np.array(th, float) / 100.0 * bdi     # gravimetric % -> vol.
        if not np.all(np.isfinite(theta)) or theta.max() > MAX_THETA \
                or theta.min() <= 0:
            skip["theta"] += 1
            continue

        if not np.isfinite(por) or not (POROSITY_RANGE[0] < por
                                        < POROSITY_RANGE[1]):
            skip["porosity"] += 1
            continue
        # Porosity computed from oven-dry bulk density sometimes lands below
        # the wettest measured water content -- common in shrink-swell clays,
        # where the dried clod has already shrunk. Saturation cannot be below
        # a measured value, so floor the anchor just above it rather than
        # discarding the layer (dropping these cost most of the clay).
        anchor = max(por * SATURATION_FRACTION, theta.max() * 1.05)
        if anchor > MAX_THETA:
            skip["porosity"] += 1
            continue
        n_measured = len(h)
        h = np.concatenate([[0.0], h])
        theta = np.concatenate([[anchor], theta])

        sand, silt, clay = row.sand_total, row.silt_total, row.clay_total
        tot = sum(x for x in (sand, silt, clay) if x is not None
                  and np.isfinite(x)) if None not in (sand, silt, clay) else 0
        if not tot or not (95 <= tot <= 105):
            skip["texture"] += 1
            continue
        sand, silt, clay = 100 * sand / tot, 100 * silt / tot, 100 * clay / tot

        try:
            popt, _ = st.fit_vg(h, theta)
        except Exception:
            skip["fit"] += 1
            continue
        tr, ts, la, ln1 = popt
        alpha, n = 10.0 ** la, 1.0 + 10.0 ** ln1
        rmse = float(np.sqrt(np.mean(
            (st.vg_theta(h, tr, ts, alpha, n) - theta) ** 2)))
        if not np.isfinite(rmse) or rmse > MAX_RMSE or n <= 1.0:
            skip["rmse"] += 1
            continue

        depth = (np.nan if row.hzn_top is None or row.hzn_bot is None
                 else (float(row.hzn_top) + float(row.hzn_bot)) / 2.0)
        rows.append(dict(
            layer_id=f"KSSL_{row.layer_key}",
            profile_id=f"KSSL_{row.pedon_key}",
            texture_class=st.usda_class(sand, silt, clay),
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=sand, silt=silt, clay=clay,
            ksat_cmh=np.nan,          # KSSL distributes none -- see docstring
            depth_cm=depth, lat=row.lat, lon=row.lon,
            rmse=rmse, n_points=n_measured, porosity=por,
            texture_lab=row.texture_lab))

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"skipped: {skip}")
    print(f"\nKSSL layers usable: {len(out)}  -> {OUT}")
    print(f"  bulk density from 1/3-bar for {n_bd33} layers, oven-dry "
          f"otherwise")
    print(f"  median retention points per curve: {out.n_points.median():.0f}")
    print(f"  median fit RMSE: {out.rmse.median():.4f}")
    print(f"  with measured Ksat: {int(out.ksat_cmh.notna().sum())}")
    print(f"  with depth: {int(out.depth_cm.notna().sum())}   "
          f"with lat/lon: {int(out.lat.notna().sum())}")
    print("\nclass distribution:")
    print(out.texture_class.value_counts().to_string())

    print(f"\nthetas > 0.75 (unphysical wet-end extrapolation): "
          f"{(out.thetas > 0.75).mean()*100:.1f}% of layers")

    # ---- unit check 1: PSDA-derived class vs KSSL's own texture_lab --------
    lab = out.copy()
    lab["lab_class"] = [lab_code_to_class(c) for c in lab.texture_lab]
    lab = lab[lab.lab_class.notna()]
    if len(lab):
        agree = (lab.lab_class == lab.texture_class)
        print(f"\nunit check 1 -- PSDA-derived class matches KSSL's own "
              f"texture_lab for {agree.sum()}/{len(lab)} layers "
              f"({agree.mean()*100:.1f}%)")

    # ---- unit check 2: fitted alpha on the same scale as GSHP -------------
    gshp = pd.read_csv(st.REFERENCE_CSV)
    print("\nunit check 2 -- class-median alpha [kPa^-1] / n, KSSL vs GSHP:")
    for cls in ["sandy loam", "loam", "silt loam", "clay loam", "clay"]:
        a, b = out[out.texture_class == cls], gshp[gshp.texture_class == cls]
        if len(a) and len(b):
            print(f"  {cls:<12s} alpha {a.alpha_kpa.median():7.3f} vs "
                  f"{b.alpha_kpa.median():7.3f}    n {a.n.median():5.2f} vs "
                  f"{b.n.median():5.2f}    ({len(a)} vs {len(b)} soils)")


if __name__ == "__main__":
    main()
