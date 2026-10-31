"""One-time preprocessing of Ksat-SSIR-DB, Brazil (Ottoni et al. 2025).

Data: Ottoni, M. V. et al. (2025). Ksat-SSIR-DB: a database of saturated
hydraulic conductivity and steady-state infiltration rate of Brazilian soils.
Revista Brasileira de Ciencia do Solo 49, e0240003,
doi:10.36783/18069657rbcs20240003. File ksat-ssir-db-2024_v3.zip (1 July
2024) from the Servico Geologico do Brasil page
https://www.sgb.gov.br/ksat-ssir-db-base-de-dados-de-condutividade-hidraulica-saturada-e-de-taxa-de-infiltracao-basica-em-solos-brasileiros
(journal article CC BY; the data licence is not stated).

Reads data/ksatssir_raw/x/KSAT-SSIR-DB-2024.xlsx (unzipped) and writes
data/ksatssir_reference.csv. Under test as an independent set
(verify_external.py ksatssir): it is not in any reference set the tool uses.

What is used, per the file's metadata sheet:
  * Rows with water content at all five suctions the database carries: 0, 6,
    10, 33 and 1500 kPa (215 of 2,579 rows). Five points is exactly
    MIN_POINTS. Gravimetric values were converted to volume by the compilers.
  * Three source theses lie at the same places as HYBRAS v1 layers already in
    the reference (all their rows within 2 km of one): Leal (2011), Rojas
    (1998) and Souza and Souza (2001). They are most likely the same samples,
    so they are left out (100 rows).
  * Texture in %, Brazilian (Embrapa) limits: sand 0.05-2 mm, as USDA.
  * Ks: the laboratory value (Kslab, mm/h) on undisturbed rings only, constant
    or falling head. The field value (Kfs) is not used.
  * Sample type from WR_sample_TP ("indeformada" = undisturbed); blank is
    "unknown". The compilers' retention-consistency note is kept in wr_note
    (for example, saturated water content far from the porosity).
"""

import os

import numpy as np
import pandas as pd

import free_m
import swcc_texture as st

RAW = os.path.join("data", "ksatssir_raw", "x", "KSAT-SSIR-DB-2024.xlsx")
OUT = os.path.join("data", "ksatssir_reference.csv")
SHEET = "KSAT-SSIR-DB - Sample Data"
MIN_POINTS = 5
MAX_RMSE = 0.03
H_KPA = {"TH0 (cm3/cm3)": 0.0, "TH6 (cm3/cm3)": 6.0, "TH10 (cm3/cm3)": 10.0,
         "TH33 (cm3/cm3)": 33.0, "TH1500 (cm3/cm3)": 1500.0}
HYBRAS_SITES = {"Leal (2011)", "Rojas (1998)", "Souza and Souza (2001)"}
NUM = ["LatitudeOR", "LongitudeOR", "Upper_Hor", "Lower_Hor", "clay (%)",
       "silt (%)", "total sand (%)", "bulk_den (g/cm3)", "porosity (cm3/cm3)",
       "org_carb (g/kg)", "Kslab (mm/h)", *H_KPA]


def read():
    """Rows with all five retention points."""
    d = pd.read_excel(RAW, SHEET)
    for c in NUM:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[d[list(H_KPA)].notna().all(axis=1)].copy()
    d["layer_id"] = "KSB_" + d.Code.astype(str)
    return d.reset_index(drop=True)


def points(r):
    return np.array(list(H_KPA.values())), np.array([r[c] for c in H_KPA])


def measured_points():
    """{layer_id: (h kPa, theta)} for every row with a full curve."""
    return {r.layer_id: points(r) for _, r in read().iterrows()}


def main():
    if not os.path.exists(RAW):
        raise SystemExit(f"{RAW} not found -- download ksat-ssir-db-2024_v3.zip "
                         "from the SGB page and unzip it into data/ksatssir_raw/x/")
    d = read()
    rows = []
    skip = dict(hybras=0, texture=0, fit=0, rmse=0)
    for _, r in d.iterrows():
        if r.Reference_simplif in HYBRAS_SITES:
            skip["hybras"] += 1
            continue
        tot = r["total sand (%)"] + r["silt (%)"] + r["clay (%)"]
        if not np.isfinite(tot) or not (95 <= tot <= 105):
            skip["texture"] += 1
            continue
        sand, silt, clay = (100 * r[c] / tot for c in
                            ("total sand (%)", "silt (%)", "clay (%)"))
        h, theta = points(r)
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
        ks = r["Kslab (mm/h)"]
        ks_ok = r.Kslab_Sample_Type == "indeformada" and np.isfinite(ks) and ks > 0
        undisturbed = r.WR_sample_TP == "indeformada"
        rows.append(dict(
            layer_id=r.layer_id,
            # No profile code in the file: one site and source is one profile.
            profile_id=f"KSB_{r.Reference_simplif}_{r.LatitudeOR:.3f}_{r.LongitudeOR:.3f}"
                       .replace(" ", ""),
            texture_class=st.usda_class(sand, silt, clay),
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=sand, silt=silt, clay=clay,
            ksat_cmh=ks / 10.0 if ks_ok else np.nan,
            depth_cm=0.5 * (r.Upper_Hor + r.Lower_Hor),
            lat=r.LatitudeOR, lon=r.LongitudeOR, state=r.State,
            # "Rodrigues et al (1991)" is also spelled with a dot.
            source_db="KsatSSIR_" + r.Reference_simplif.replace(" ", "").replace(".", ""),
            oc=r["org_carb (g/kg)"] / 10.0, bd=r["bulk_den (g/cm3)"],
            porosity=r["porosity (cm3/cm3)"], texture_lab=r.Texture_en,
            soil_class=r.Soil_Class_u2NC,
            wr_note=r.WRetention_Consistency if isinstance(
                r.WRetention_Consistency, str) else "",
            rmse=rmse, n_points=len(h), **fm,
            sample_type="undisturbed" if undisturbed else "unknown",
            sample_type_source="Ksat-SSIR-DB WR_sample_TP"))
    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"rows with five points: {len(d)}; wrote {OUT}: {len(out)} rows, "
          f"{out.profile_id.nunique()} profiles, {out.source_db.nunique()} sources "
          f"(skipped: {skip})")
    print(f"  by state: {out.state.value_counts().to_dict()}")
    print(f"  Ks: {out.ksat_cmh.notna().sum()} (median {out.ksat_cmh.median():.1f} cm/h); "
          f"undisturbed: {(out.sample_type == 'undisturbed').sum()}; "
          f"median fit RMSE {out.rmse.median():.4f}")
    lab = out.texture_lab.astype(str).str.lower()
    print(f"  class from fractions matches the file's class: "
          f"{(lab == out.texture_class).sum()}/{len(out)}")
    print(out.texture_class.value_counts().to_string())


if __name__ == "__main__":
    main()
