"""One-time preprocessing of the Laikipia, Kenya soils (Willard et al.).

Willard, L. L., Muñoz-Carpena, R., Venort, T., Gitonga, J., Maltais-Landry,
G., Caylor, K. K. and Palm, C. A. A high-resolution comprehensive database to
evaluate agricultural impacts at edge of savanna in Kenya. Under review,
Scientific Data.

Reads data/willard_raw/SoilsData.csv (the revised long-format release of
2026-09-18, one row per sample x measurement) and, for the sandbox segment
only, data/willard_raw/SoilsData_ExactLocation.csv (the original laboratory
records), and data/willard_raw/corrected/
SWCC_shrinkage_summary_revised_2026-10-16.csv (the authors' revised curves,
see "Final curves" below). Writes data/willard_reference.csv. The data are
the project team's own: the derived table is distributed (2026-10-31); the raw
files stay local, since the original records carry exact site locations.

What is used, per the data's ReadMe:
  * 86 point samples: grid locations 1 and 4 in each of 43 fields
    (smallholder, large farm or adjacent bush) in Laikipia, Nyeri and Meru,
    January-February 2022; undisturbed cores at 15 cm.
  * Retention: 13 points, pF 0 to 4.2. The sandbox (pF 0-2) did not
    saturate the cores, and the released file raises those five points by
    one constant per sample so that pF 0 sits at 95 % of porosity, leaving
    pF 2.15-4.2 as measured. That opens a step of median 0.21 between pF 2
    and pF 2.15 -- half of each curve's whole water loss -- which no clay
    shows.

    The ORIGINAL laboratory values (SoilsData_ExactLocation.csv) are used
    here for pF 0-2 instead: unshifted, they run continuously into the
    sand-kaolinite segment (pF 2 at 27.4 % against pF 2.15 at 26.5 %). The
    released file supplies pF 2.15-4.2, which carries the eight drying-
    sequence corrections the authors made, and the revised bulk density and
    Ks.

    The whole curve is then multiplied by one factor per sample,
    0.95 x porosity / theta(pF 0), about 1.5. Rationale (2026-09-20): the
    original curves have the right shape but the wrong level -- as measured
    the tool reads them as loams (2 % exact), and the field observation is
    that these soils are heavy. Corrections confined to the sandbox, flat
    or tapered, leave the middle and dry range too dry to be clay and reach
    only ~53 % exact; the pressure-plate values are low against reference
    clays as well (0.186 against 0.284 at wilting point). One factor per
    sample keeps the measured shape, closes the junction step to 0.014, and
    puts the dry end at 0.283 -- reference clays sit at 0.284, which nothing
    in the procedure aimed at. It reads as a volumetric-basis error (ring
    volume or the bulk density used in the conversion) rather than an
    apparatus fault; true shrinkage would need a factor that grows with
    suction, and would also put these curves on a different basis from the
    core-volume theta of every source they are compared against.

    Shrinking density (2026-10-12, project decision): these are vertic
    clays, which shrink as they dry, so the curve -- volumetric at one wet
    density per sample, the BD1 = 2.65 x (1 - porosity) behind the scaling
    above -- is converted point by point to a density that follows the
    water content, as for the Tomejil vertisol (prepare_andalusia.py):
    w = theta / BD1, BD(w) = min(1.75, 1 / (1/2.65 + w)) (normal shrinkage
    up to the median oven-dry clod density of KSSL Vertisols like these),
    theta = w * BD(w). The wet end barely moves (x1.03 at pF 0) and the dry
    end rises by about half (x1.47 at pF 4.2). As a new source the tool then
    reads 97 % of the clays as clay with 51.7 % clay predicted (measured
    51.6 %), against 81 % and 43 % with one density (Vertisol_shrinkage/
    test_caps.py). This is the "factor that grows with suction" weighed
    above: a model of shrinkage, not a measured curve.

    Final curves (2026-10-16, project decision): the reference takes the
    authors' own revision instead, the column wet_mass_1.6_shrinkage of
    SWCC_shrinkage_summary_revised_2026-10-16.csv (L. Willard), which the
    Scientific Data paper will publish: theta from the laboratory's wet
    masses, converted with the same normal-shrinkage density capped at
    1.6 g/cm3. The 16 samples the authors flag with dispersion = 1 (their
    hydrometer clay is underestimated, from volcanic characteristics) are
    left out, so 70 samples remain. As a new source the tool reads 93 % of
    them to the exact class and 98.5 % of the clays as clay (cap 1.75 above:
    89 % and 94 %; their wet_mass_shrink_group95 column: 47 % and 50 %).
    The scaled curve above (scaled_curve) is kept for the Vertisol_shrinkage
    memo.

    Samples whose bulk density is missing take their soil cluster's median
    porosity. h_kPa = 10**pF * 0.0980665, pF 0 taken as h = 0; theta is
    volume %. Blank points are dropped. The samples were discarded, so the
    sandbox cannot be rerun.
  * Texture: sieve plus hydrometer, USDA limits (sand > 50 um). The class is
    re-derived from the fractions, as for every other source.
  * Ks: the Saturo dual-head infiltrometer value is the sample's Ks
    (ksat_cmh; project decision, 2026-10-16, 26 of the 70 samples). Against
    the final curves its physical estimate sits within a factor of 10 for
    92 % of samples, the kNN for 54 %. Until then no Ks entered the
    reference (decision of 2026-09-18). All three conductivities were
    measured in the field in the dry season, when these vertic clays were at or below wilting point and
    cracked: the permeameter's median is 178 cm/h on the drier clay-rich
    samples against 10-14 cm/h elsewhere, and every method rises with clay.
    The reference's Ks is mostly saturated, swollen lab cores. The field
    values are kept, m/s -> cm/h, as ksat_mpd_cmh (modified Philip-Dunne
    permeameter; the value the authors reject, 3LN06324, is left out),
    ksat_saturo_cmh (dual-head infiltrometer) and k_minidisk_cmh
    (unsaturated, 2 cm suction).
  * Bulk density: mean of the KALRO and UF cores of the sample, the value
    the authors use for porosity; a core they reject is left out. Kept as
    bd_core. The bd column, which the tool matches on when a user supplies a
    bulk density, is the density at 33 kPa on the sample's shrinkage curve
    (min(1.6, 2.65 x (1 - theta33)): on the normal-shrinkage line the soil
    stays saturated, so theta = 1 - BD/2.65 until the cap; median ~1.4), as for Tomejil and Lebrija and as KSSL reports vertisols:
    next to the shrinking-density curves the wet core density (median 1.06)
    misled the classifier (new source with depth and bulk density: 53.5 %
    exact, 62 % of clays read as clay, against 80 % and 93 % curve only).
  * Organic carbon: total carbon of each field's composite (location 0),
    carried to both samples of the field; the soils are not calcareous.
    Organic matter = 1.724 x organic carbon.
  * Coordinates: rounded by the authors to 0.01 degree; none for large farms.
  * Sample ids: 33 samples of fields LT077 and SN078 carry soil cluster 2 in
    the original file and 3 in the release; the release is followed.
"""

import os

import numpy as np
import pandas as pd

import free_m
import swcc_texture as st

RAW = os.path.join("data", "willard_raw", "SoilsData.csv")
ORIG = os.path.join("data", "willard_raw", "SoilsData_ExactLocation.csv")
OUT = os.path.join("data", "willard_reference.csv")
CM_TO_KPA = 0.0980665
M_S_TO_CM_H = 100.0 * 3600.0
DEPTH_CM = 15.0
OM_PER_OC = 1.724
MIN_POINTS = 5
MAX_RMSE = 0.03
REJECTED = 20            # QAQC_stat_flag: do not use
PF_SANDBOX = (0.0, 2.0)  # the sandbox segment, taken from ORIG and rescaled
PARTICLE_DENSITY = 2.65
SATURATION_FRACTION = 0.95   # theta(pF 0) as a fraction of porosity
SHRINK_CAP = 1.75            # g/cm3, cap of scaled_curve (see docstring)
REVISED = os.path.join("data", "willard_raw", "corrected",
                       "SWCC_shrinkage_summary_revised_2026-10-16.csv")
REVISED_COL = "wet_mass_1.6_shrinkage"   # the authors' final curve, volume %
REVISED_CAP = 1.6            # g/cm3, its shrinkage cap


def read():
    d = pd.read_csv(RAW, low_memory=False)
    d["Text"] = d.Value.astype(str).str.strip()      # class labels
    d["Value"] = pd.to_numeric(d.Value, errors="coerce")
    return d[d.QAQC_stat_flag != REJECTED]


def original_sandbox():
    """{SampleID: {pF: theta fraction}} for pF 0-2, from the original records.

    Keyed by the release's sample id: 33 samples of LT077 and SN078 differ
    from the original only in the leading soil-cluster digit.
    """
    o = pd.read_csv(ORIG, low_memory=False).set_index("SampleID")
    cols = {c: float(c[2:].split("_")[0]) for c in o.columns
            if c.startswith("pF") and c.endswith("_percent_soilmoisture")}
    out = {}
    for sid, row in o.iterrows():
        vals = {pf: pd.to_numeric(row[c], errors="coerce") / 100.0
                for c, pf in cols.items() if pf <= PF_SANDBOX[1]}
        for lead in "12345":                       # the cluster-digit shift
            out.setdefault(lead + str(sid)[1:], vals)
        out[sid] = vals
    return out


def porosity(d):
    """Porosity per sample from the mean of the two bulk densities, with the
    soil cluster's median where a sample has none."""
    w = d[d.SampleType == "point_sample"].pivot_table(
        index="SampleID", columns="Measurement", values="Value",
        aggfunc="first")
    bd = pd.concat([w.get("BulkDensity_KALRO"), w.get("BulkDensity_UF")],
                   axis=1).mean(axis=1)
    por = 1.0 - bd / PARTICLE_DENSITY
    cluster = pd.Series([str(i)[0] for i in por.index], index=por.index)
    return por.fillna(cluster.map(por.groupby(cluster).median()))


def curve(d=None):
    """The authors' final curves: {SampleID: (h kPa, theta, BD at 33 kPa)},
    samples flagged dispersion = 1 left out."""
    r = pd.read_csv(REVISED)
    bad = set(r.loc[r.dispersion == 1, "SampleID"])
    out = {}
    for s, g in r.dropna(subset=[REVISED_COL]).groupby("SampleID"):
        if s in bad:
            continue
        g = g.sort_values("pF")
        h = np.where(g.pF == 0, 0.0, 10.0 ** g.pF * CM_TO_KPA)
        th = g[REVISED_COL].to_numpy() / 100.0
        pos = h > 0
        th33 = np.interp(np.log10(33.0), np.log10(h[pos]), th[pos])
        out[s] = (h, th, float(min(REVISED_CAP, PARTICLE_DENSITY * (1.0 - th33))))
    return out


def scaled_curve(d):
    """Retention points per sample: {SampleID: (h kPa, theta fraction)}.

    pF 0-2 come from the original records, pF 2.15-4.2 from the release, the
    whole curve is scaled so that pF 0 sits at 95 % of porosity, and it is
    then converted to a shrinking density (see the module docstring).
    """
    r = d[d.Measurement.str.match(r"^pF[\d.]+_soilmoisture$")].dropna(
        subset=["Value"])
    pf = r.Measurement.str.extract(r"^pF([\d.]+)_")[0].astype(float)
    r = r.assign(pf=pf, theta=r.Value / 100.0)
    sandbox, por = original_sandbox(), porosity(d)
    out = {}
    for s, g in r.groupby("SampleID"):
        orig = sandbox.get(s)
        if orig is None:
            continue
        theta = {p: (orig[p] if p <= PF_SANDBOX[1] and p in orig else t)
                 for p, t in zip(g.pf, g.theta)}
        theta = {p: t for p, t in theta.items() if np.isfinite(t) and t > 0}
        p0 = theta.get(0.0)
        target = SATURATION_FRACTION * por.get(s, np.nan)
        if not p0 or not np.isfinite(target):
            continue
        f = target / p0
        pfs = np.array(sorted(theta))
        h = np.where(pfs == 0, 0.0, 10.0 ** pfs * CM_TO_KPA)
        th1 = np.array([theta[p] for p in pfs]) * f        # one density, BD1
        w = th1 / (PARTICLE_DENSITY * (1.0 - por[s]))       # gravimetric
        bd = np.minimum(SHRINK_CAP, 1.0 / (1.0 / PARTICLE_DENSITY + w))
        # density at 33 kPa on the same shrinkage curve (w interpolated in
        # log suction), the bulk density the reference stores for vertisols
        pos = h > 0
        w33 = np.interp(np.log10(33.0), np.log10(h[pos]), w[pos])
        bd33 = float(min(SHRINK_CAP, 1.0 / (1.0 / PARTICLE_DENSITY + w33)))
        out[s] = (h, w * bd, bd33)
    return out


def measured_points():
    """{layer_id: (h kPa, theta)} for every sample, for tests that refit."""
    return {f"WL_{s}": (h, th) for s, (h, th, _) in curve(read()).items()}


def main():
    if not os.path.exists(RAW):
        raise SystemExit(f"{RAW} not found")
    d = read()
    pts = curve(d)
    wide = d[d.SampleType == "point_sample"].pivot_table(
        index="SampleID", columns="Measurement", values="Value", aggfunc="first")
    site = d.drop_duplicates("SampleID").set_index("SampleID")
    oc = (d[(d.Measurement == "TotalCarbon")]
          .drop_duplicates("SiteID").set_index("SiteID").Value)

    rows = []
    skip = dict(points=0, texture=0, fit=0, rmse=0)
    for s, (h, theta, bd33) in pts.items():
        if len(h) < MIN_POINTS:
            skip["points"] += 1
            continue
        w = wide.loc[s]
        sand, silt, clay = w.Sand, w.Silt, w.Clay
        tot = sand + silt + clay
        if not np.isfinite(tot) or not (95 <= tot <= 105):
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
        fm = free_m.fit_free_m(h, theta,
                               fallback=free_m.mualem_fallback(tr, ts, alpha, n))

        info = site.loc[s]
        c = oc.get(info.SiteID, np.nan)
        bd = np.nanmean([w.get("BulkDensity_KALRO", np.nan),
                         w.get("BulkDensity_UF", np.nan)])
        rows.append(dict(
            layer_id=f"WL_{s}",
            # The two samples of a field share a site and a depth, so they
            # are grouped to keep them out of each other's hold-out folds.
            profile_id=f"WL_{info.SiteID}",
            texture_class=st.usda_class(sand, silt, clay),
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=sand, silt=silt, clay=clay,
            ksat_cmh=w.get("Ksat_Saturo", np.nan) * M_S_TO_CM_H,
            ksat_mpd_cmh=w.get("Ksat_MPD", np.nan) * M_S_TO_CM_H,
            ksat_saturo_cmh=w.get("Ksat_Saturo", np.nan) * M_S_TO_CM_H,
            k_minidisk_cmh=w.get("K_MiniDisk", np.nan) * M_S_TO_CM_H,
            depth_cm=DEPTH_CM, lat=info.Latitude, lon=info.Longitude,
            source_db="Willard_Laikipia", oc=c, om=OM_PER_OC * c,
            porosity=w.get("Porosity", np.nan), bd=bd33, bd_core=bd,
            rmse=rmse, n_points=len(h), **fm,
            sample_type="undisturbed",
            sample_type_source="Willard et al.: retention on undisturbed "
                               "cores (KALRO); Ks in situ (Saturo)"))

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"samples with a curve: {len(pts)}")
    print(f"reference layers: {len(out)}   (skipped: {skip})")
    print(f"  fields:         {out.profile_id.nunique()}  "
          f"(with coordinates: {out.lat.notna().sum()} samples)")
    print(f"  Ks (Saturo): {out.ksat_cmh.notna().sum()} (median {out.ksat_cmh.median():.1f} cm/h)")
    print(f"  field Ks kept aside: permeameter {out.ksat_mpd_cmh.notna().sum()} "
          f"(median {out.ksat_mpd_cmh.median():.1f} cm/h), Saturo "
          f"{out.ksat_saturo_cmh.notna().sum()}, mini-disk "
          f"{out.k_minidisk_cmh.notna().sum()}")
    print(f"  with bd:        {out.bd.notna().sum()} (median at 33 kPa {out.bd.median():.2f}; "
          f"core {out.bd_core.median():.2f})")
    print(f"  with oc / om:   {out.oc.notna().sum()}  (median om "
          f"{out.om.median():.1f} %)")
    print(f"  median fit RMSE: {out.rmse.median():.4f}")
    lab = (d[d.Measurement == "TextureClass"].drop_duplicates("SampleID")
           .set_index("SampleID").Text.str.lower())
    agree = out.texture_class.to_numpy() == lab.reindex(
        [i[3:] for i in out.layer_id]).to_numpy()
    print(f"  class from fractions matches the published label: "
          f"{agree.sum()}/{len(out)}")
    print(out.texture_class.value_counts().to_string())


if __name__ == "__main__":
    main()
