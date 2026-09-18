"""One-time preprocessing of the GSHP database.

Reads the raw GSHP file (data/WRC_dataset_surya_et_al_2021_final.csv,
Gupta et al. 2022, https://doi.org/10.5281/zenodo.6640246) and writes a
compact per-layer reference table data/gshp_reference.csv used by
swcc_texture.py.

Units in the output table:
    alpha_kpa  [kPa^-1]   (GSHP alpha is in m^-1 of water head; 1 m = 9.80665 kPa)
    thetar, thetas [m3/m3]
    sand, silt, clay [%]  (renormalized to sum to 100)
    ksat_cmh   [cm/h]     (GSHP ksat_lab / ksat_field are in cm/day, except
                           Florida_database, whose values are already cm/h)
"""

import pandas as pd

import free_m

M_HEAD_TO_KPA = 9.80665
# lab_head_m carries a 2.48e29 sentinel for missing values.
MAX_HEAD_M = 1e5

USDA_CLASSES = [
    "sand", "loamy sand", "sandy loam", "loam", "silt", "silt loam",
    "sandy clay loam", "clay loam", "silty clay loam", "sandy clay",
    "silty clay", "clay",
]


def main():
    df = pd.read_csv("data/WRC_dataset_surya_et_al_2021_final.csv",
                     low_memory=False, encoding="latin-1")
    lay = df.drop_duplicates("layer_id").copy()

    for c in ["alpha", "n", "thetar", "thetas", "sand_tot_psa", "silt_tot_psa",
              "clay_tot_psa", "ksat_lab", "ksat_field", "hzn_top", "hzn_bot",
              "latitude_decimal_degrees", "longitude_decimal_degrees"]:
        lay[c] = pd.to_numeric(lay[c], errors="coerce")

    lay["tex_psda"] = lay["tex_psda"].astype(str).str.strip().str.lower()
    lay = lay[lay["tex_psda"].isin(USDA_CLASSES)]
    lay = lay[lay["data_flag"] == "good quality estimate"]
    lay = lay[(lay["n"] > 1.0) & (lay["alpha"] > 0)]

    out = pd.DataFrame({
        "layer_id": lay["layer_id"],
        "profile_id": lay["profile_id"],
        "texture_class": lay["tex_psda"],
        "alpha_kpa": lay["alpha"] / M_HEAD_TO_KPA,
        "n": lay["n"],
        "thetar": lay["thetar"],
        "thetas": lay["thetas"],
    })

    # Particle fractions: keep only triples that roughly close to 100 %,
    # then renormalize exactly.
    s = lay["sand_tot_psa"] + lay["silt_tot_psa"] + lay["clay_tot_psa"]
    ok = s.between(95, 105)
    for col, raw in [("sand", "sand_tot_psa"), ("silt", "silt_tot_psa"),
                     ("clay", "clay_tot_psa")]:
        out[col] = (lay[raw] / s * 100).where(ok)

    # Ksat: prefer lab, fall back to field; convert cm/day -> cm/h. Florida is
    # the exception: its values are already cm/h. The Florida Soil
    # Characterization metadata gives KSat in cm/hr, and GSHP's Florida medians
    # read as cm/h match the class means (sand 24, sandy loam 1.1, clay 0.14;
    # Carsel & Parrish 29.7, 4.4, 0.20), while read as cm/day they fall 24x low.
    ksat = lay["ksat_lab"].fillna(lay["ksat_field"])
    per_day = pd.Series(24.0, index=lay.index).where(
        lay["source_db"] != "Florida_database", 1.0)
    out["ksat_cmh"] = ksat.where(ksat > 0) / per_day

    # Sample mid-depth (cm), used as an optional covariate.
    out["depth_cm"] = (lay["hzn_top"] + lay["hzn_bot"]) / 2.0

    # Oven-dry bulk density (g/cm3), an optional user covariate.
    bd = pd.to_numeric(lay["db_od"], errors="coerce")
    out["bd"] = bd.where((bd > 0.1) & (bd < 2.3))

    # Coordinates, used to test how much region-matched reference data matters
    # (see verify_region.py).
    out["lat"] = lay["latitude_decimal_degrees"]
    out["lon"] = lay["longitude_decimal_degrees"]

    # Contributing database. GSHP is dominated by one source (Florida is ~57 %
    # of the curated layers), so profile-level leave-one-out still lets a
    # target be matched against siblings measured by the same lab with the same
    # protocol. verify_source_blocked.py holds out a whole source to measure
    # how much that inflates the accuracy estimates.
    out["source_db"] = lay["source_db"]

    # Standard errors of GSHP's own vG fit, as RELATIVE errors: alpha spans
    # orders of magnitude, so absolute se scales with alpha and is not
    # comparable across soils. Used to down-weight poorly identified reference
    # layers (see verify_gbm.py). Median rse_alpha is 0.11 and only 0.9 % of
    # layers exceed 1.0, so the reference is mostly well constrained.
    out["rse_alpha"] = pd.to_numeric(lay["se_alpha"], errors="coerce") \
        / lay["alpha"]
    out["rse_n"] = pd.to_numeric(lay["se_n"], errors="coerce") \
        / (lay["n"] - 1.0)

    # Organic carbon (%), an optional covariate. Coverage is only ~15 % in
    # GSHP, against ~85 % in KSSL -- see verify_om.py for whether it earns the
    # reference rows it costs.
    oc = pd.to_numeric(lay["oc"], errors="coerce")
    out["oc"] = oc.where((oc >= 0) & (oc < 60))

    # Sample type: how the WET end of the retention curve was measured -- an
    # undisturbed sample (core, clod or in situ) versus repacked or sieved
    # material. A dry end measured on sieved soil is standard practice for
    # undisturbed samples and does not change the label. GSHP's own
    # disturbed_undisturbed field is used, except for Florida_database, which
    # GSHP leaves "unknown": one of the dataset's authors (W. G. Harris, pers.
    # comm., 2026-09-12) recalls the rings for conductivity and water release
    # being collected in the field, per horizon, as undisturbed samples, with
    # no repacking, and the release curves measured at many pressure steps.
    stype = lay["disturbed_undisturbed"].fillna("unknown").astype(str) \
        .str.strip().str.lower()
    fl = lay["source_db"] == "Florida_database"
    out["sample_type"] = stype.where(~fl, "undisturbed")
    out["sample_type_source"] = (
        pd.Series("GSHP disturbed_undisturbed field (Gupta et al. 2022)",
                  index=out.index)
        .where(stype != "unknown",
               "GSHP lists it as unknown; the source gives no preparation")
        .where(~fl, "Florida: field rings per horizon, undisturbed "
                    "(W. G. Harris, pers. comm.)"))

    # Unconstrained-m parameters. GSHP publishes only the Mualem-constrained
    # fit, but it also ships the measured (h, theta) points those were fitted
    # to, so m can be freed by refitting. Layers with too few points keep the
    # published parameters with m = 1 - 1/n (see free_m.py).
    pts = df[["layer_id", "lab_head_m", "lab_wrc"]].copy()
    for c in ("lab_head_m", "lab_wrc"):
        pts[c] = pd.to_numeric(pts[c], errors="coerce")
    pts = pts[(pts.lab_head_m >= 0) & (pts.lab_head_m < MAX_HEAD_M)
              & (pts.lab_wrc > 0) & (pts.lab_wrc <= 1)]
    by_layer = {k: v for k, v in pts.groupby("layer_id")}

    rows, n_free = [], 0
    for r in out.itertuples():
        fb = free_m.mualem_fallback(r.thetar, r.thetas, r.alpha_kpa, r.n)
        g = by_layer.get(r.layer_id)
        res = fb
        if g is not None:
            res = free_m.fit_free_m(g.lab_head_m.to_numpy() * M_HEAD_TO_KPA,
                                    g.lab_wrc.to_numpy(), fallback=fb)
        n_free += res is not fb
        rows.append(res)
    for c in free_m.FREE_COLS:
        out[c] = [x[c] for x in rows]
    print(f"  free-m refit:   {n_free} of {len(out)} layers "
          f"({n_free / len(out) * 100:.1f} %); rest keep m = 1 - 1/n")

    out.to_csv("data/gshp_reference.csv", index=False)
    print(f"reference layers: {len(out)}")
    print(f"  with fractions:  {out['sand'].notna().sum()}")
    print(f"  with ksat:       {out['ksat_cmh'].notna().sum()}")
    print(f"  with depth:      {out['depth_cm'].notna().sum()}")
    print(f"  sample type:     {out.sample_type.value_counts().to_dict()}")
    print(out["texture_class"].value_counts())


if __name__ == "__main__":
    main()
