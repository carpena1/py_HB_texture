"""One-time preprocessing of the Kansas Mesonet soils (Parker et al. 2022).

Data: Parker, N., Kluitenberg, G. J., Redmond, C. and Patrignani, A. (2022).
A database of soil physical properties for the Kansas Mesonet. Soil Science
Society of America Journal 86(6), 1495-1508, doi:10.1002/saj2.20465. The 2026
summary file (QA/QC corrections, qaqc_log tab) is public at
https://soilwater.github.io/kansas-soil-props/ (repository soilwater/
kansas-soil-props, MIT licence); A. Patrignani pointed us to it (2026-10-30).

Reads data/kansas_raw/kansas_mesonet_soil_database_summary_2026.xlsx and
writes data/kansas_reference.csv. Under test as an independent set
(verify_external.py kansas): it is not in any reference set the tool uses.

What is used, per the file's metadata tab and Parker et al. (2022):
  * 320 undisturbed cores (Eijkelkamp ring sampler) from 40 Mesonet stations,
    two cores at each of 5, 10, 20 and 50 cm; 4 cores could not be taken.
  * Retention on the same core: saturation (several days in CaCl2 solution),
    2 and 5 kPa on a tension table, 10, 33 and 70 kPa in Tempe cells, 1500 kPa
    on a pressure plate, water by oven drying. The saturation point is dropped
    where the qaqc_log flags it: above the porosity by more than 0.02, the core
    swelled after wetting (13 cores, most at Clay Center).
  * Texture by hydrometer, USDA limits, on ONE core of each station and depth
    (160 cores); the file gives the other core the same class but no
    fractions. Only cores with measured fractions enter the table; the class
    is re-derived from them. measured_points() returns every core, for the Ks
    tests that need no texture.
  * Ks: closed-path permeameter, constant head, on the core. The lab notes
    many cores: impermeable or never saturated (Ks 0), gaps between soil and
    ring, partly filled rings, and flow along roots, cracks or macropores.
    ksat_cmh keeps only cores with a positive Ks and no note; the note and
    the raw value stay in ksat_note and ksat_all_cmh.
  * Bulk density: dry mass over ring volume; porosity from it and the
    pycnometer particle density.
  * Organic matter by loss on ignition; organic carbon = OM / 1.724. (Total
    carbon by combustion includes carbonate in the western soils, so it is not
    used.)
"""

import os

import numpy as np
import pandas as pd

import free_m
import swcc_texture as st

RAW = os.path.join("data", "kansas_raw",
                   "kansas_mesonet_soil_database_summary_2026.xlsx")
OUT = os.path.join("data", "kansas_reference.csv")
OM_PER_OC = 1.724
MIN_POINTS = 5
MAX_RMSE = 0.03
SAT_EXCESS = 0.02        # qaqc_log: theta_sat above porosity by more than this
H_KPA = {"vwc_2kPa": 2.0, "vwc_5kPa": 5.0, "vwc_10kPa": 10.0,
         "vwc_33kPa": 33.0, "vwc_70kPa": 70.0, "vwc_1500kPa": 1500.0}
KEY = ["station_name", "ring_number"]


def layer_id(station, ring):
    return f"KS_{station.replace(' ', '')}_{int(ring)}"


def read():
    """One row per core: general and hydraulic properties joined."""
    g = pd.read_excel(RAW, "general_properties")
    h = pd.read_excel(RAW, "hydraulic_properties")
    d = g.merge(h.drop(columns=["county", "latitude", "longitude", "core_number",
                                "top_depth", "bottom_depth", "nominal_depth"]),
                on=KEY, how="inner", validate="one_to_one")
    d = d[d[list(H_KPA)].notna().all(axis=1)].reset_index(drop=True)
    d["layer_id"] = [layer_id(s, r) for s, r in zip(d.station_name, d.ring_number)]
    return d


def points(r):
    h, th = list(H_KPA.values()), [r[c] for c in H_KPA]
    if np.isfinite(r.vwc_sat) and not r.vwc_sat > r.porosity + SAT_EXCESS:
        h, th = [0.0] + h, [r.vwc_sat] + th
    return np.array(h), np.array(th)


def measured_points():
    """{layer_id: (h kPa, theta)} for every core with a full curve."""
    return {r.layer_id: points(r) for _, r in read().iterrows()}


def main():
    if not os.path.exists(RAW):
        raise SystemExit(f"{RAW} not found -- download it from "
                         "https://soilwater.github.io/kansas-soil-props/")
    d = read()
    rows = []
    skip = dict(texture=0, fit=0, rmse=0)
    for _, r in d.iterrows():
        tot = r.sand + r.silt + r.clay
        if not np.isfinite(tot) or not (95 <= tot <= 105):
            skip["texture"] += 1
            continue
        sand, silt, clay = (100 * r[c] / tot for c in ("sand", "silt", "clay"))
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
        note = r.Ksat_notes if isinstance(r.Ksat_notes, str) else ""
        ks = r.Ksat_cm_per_hr
        clean = np.isfinite(ks) and ks > 0 and not note
        om = r.OM
        rows.append(dict(
            layer_id=r.layer_id,
            # The two cores of a station share its site and soil, so they are
            # grouped to keep them out of each other's hold-out folds.
            profile_id=f"KS_{r.station_name.replace(' ', '')}",
            texture_class=st.usda_class(sand, silt, clay),
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=sand, silt=silt, clay=clay,
            ksat_cmh=ks if clean else np.nan,
            ksat_all_cmh=ks, ksat_note=note,
            depth_cm=0.5 * (r.top_depth + r.bottom_depth),
            lat=r.latitude, lon=r.longitude, source_db="Kansas_Mesonet",
            oc=om / OM_PER_OC, om=om, porosity=r.porosity, bd=r.bulk_density,
            rmse=rmse, n_points=len(h), **fm,
            sample_type="undisturbed",
            sample_type_source="Parker et al. (2022): retention and Ks on the "
                               "same undisturbed core"))
    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"cores with a full curve: {len(d)}; saturation point dropped "
          f"(swelling): {int((d.vwc_sat > d.porosity + SAT_EXCESS).sum())}")
    print(f"wrote {OUT}: {len(out)} cores from {out.profile_id.nunique()} "
          f"stations (skipped: {skip})")
    print(f"  Ks kept: {out.ksat_cmh.notna().sum()} (median "
          f"{out.ksat_cmh.median():.2f} cm/h); noted or zero: "
          f"{out.ksat_cmh.isna().sum()}")
    print(f"  median fit RMSE: {out.rmse.median():.4f}")
    lab = d.set_index("layer_id").textural_class.str.lower()
    agree = out.texture_class.to_numpy() == lab.reindex(out.layer_id).to_numpy()
    print(f"  class from fractions matches the file's class: "
          f"{agree.sum()}/{len(out)}")
    print(out.texture_class.value_counts().to_string())


if __name__ == "__main__":
    main()
