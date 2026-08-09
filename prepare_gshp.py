"""One-time preprocessing of the GSHP database.

Reads the raw GSHP file (data/WRC_dataset_surya_et_al_2021_final.csv,
Gupta et al. 2022, https://doi.org/10.5281/zenodo.6640246) and writes a
compact per-layer reference table data/gshp_reference.csv used by
swcc_texture.py.

Units in the output table:
    alpha_kpa  [kPa^-1]   (GSHP alpha is in m^-1 of water head; 1 m = 9.80665 kPa)
    thetar, thetas [m3/m3]
    sand, silt, clay [%]  (renormalized to sum to 100)
    ksat_cmh   [cm/h]     (GSHP ksat_lab / ksat_field are in cm/day)
"""

import pandas as pd

M_HEAD_TO_KPA = 9.80665

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

    # Ksat: prefer lab, fall back to field; convert cm/day -> cm/h.
    ksat = lay["ksat_lab"].fillna(lay["ksat_field"])
    out["ksat_cmh"] = ksat.where(ksat > 0) / 24.0

    # Sample mid-depth (cm), used as an optional covariate.
    out["depth_cm"] = (lay["hzn_top"] + lay["hzn_bot"]) / 2.0

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

    out.to_csv("data/gshp_reference.csv", index=False)
    print(f"reference layers: {len(out)}")
    print(f"  with fractions:  {out['sand'].notna().sum()}")
    print(f"  with ksat:       {out['ksat_cmh'].notna().sum()}")
    print(f"  with depth:      {out['depth_cm'].notna().sum()}")
    print(out["texture_class"].value_counts())


if __name__ == "__main__":
    main()
