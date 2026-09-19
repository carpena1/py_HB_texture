"""One-time preprocessing of the CSIRO Boorowa Farm soil physics data.

CSIRO (2025). Characterisation of Soil Physical Properties Across Seven
Representative Sites at Boorowa Farm. CSIRO Data Access Portal,
doi:10.25919/rxts-tr35 (CC BY-NC-SA 4.0: raw and derived data stay local).

Reads data/boorowa_raw/SWRC_Analyses_BoorowaFarm_2020.xlsx and writes
data/boorowa_reference.csv (git-ignored).

  * Seven sites at the CSIRO Boorowa Agricultural Research Station, NSW,
    each a different soil type, sampled 2020.
  * Retention measured at the CSIRO Soil Physics Facility (McKenzie et al.
    2002): suction tables (method 504.01) at 10, 30, 50, 100, 340 and 600 cm,
    pressure plates (504.02) at 5 and 15 bar; volumetric water content.
    Intact cores (0-10 and 10-20 cm at every site, deeper at some) are
    undisturbed; repacked samples (to 140 cm) are disturbed.
  * Texture by hydrometer on the fine earth (Gee and Bauder 1986), per site
    and depth; organic carbon as TOC. No Ks.
  * Curves are refitted with the tool's own routine; the file's own fits are
    not used.
"""

import os

import numpy as np
import pandas as pd

import free_m
import swcc_texture as st

RAW = os.path.join("data", "boorowa_raw", "SWRC_Analyses_BoorowaFarm_2020.xlsx")
OUT = os.path.join("data", "boorowa_reference.csv")
CM_TO_KPA = 0.0980665
H_CM = np.array([10, 30, 50, 100, 340, 600, 5098.58, 15295.7])
VOL_COLS = list(range(20, 28))       # volumetric water content at H_CM
MIN_POINTS = 5
MAX_RMSE = 0.03


def depth_key(s):
    return str(s).strip().rstrip(".")


def layer_id(site, depth, condition):
    return f"BW_{site}_{depth}_{condition[0]}"


def samples():
    d = pd.read_excel(RAW, sheet_name="Data", header=None)
    d = d[pd.to_numeric(d[1], errors="coerce").notna()]
    out = []
    for _, r in d.iterrows():
        theta = pd.to_numeric(r[VOL_COLS], errors="coerce").to_numpy(float)
        ok = np.isfinite(theta)
        out.append(dict(site=int(r[1]), condition=str(r[2]).strip(),
                        depth=depth_key(r[3]), bd=pd.to_numeric(r[4], errors="coerce"),
                        h=H_CM[ok] * CM_TO_KPA, theta=theta[ok]))
    return out


def measured_points():
    return {layer_id(s["site"], s["depth"], s["condition"]): (s["h"], s["theta"])
            for s in samples()}


def main():
    if not os.path.exists(RAW):
        raise SystemExit(f"{RAW} not found")
    tex = pd.read_excel(RAW, sheet_name="soil_attributes")
    tex["depth"] = tex.depth.map(depth_key)
    tex = tex.set_index(["PAWC_site", "depth"])
    loc = pd.read_excel(RAW, sheet_name="location info").set_index("Site")
    rows = []
    skip = dict(points=0, texture=0, fit=0, rmse=0)
    for s in samples():
        h, theta = s["h"], s["theta"]
        if len(h) < MIN_POINTS:
            skip["points"] += 1
            continue
        key = (s["site"], s["depth"])
        if key not in tex.index:
            skip["texture"] += 1
            continue
        t = tex.loc[key]
        frac = np.array([t.sand_fs, t.silt_fs, t.clay_fs], float) * 100
        tot = frac.sum()
        if not np.isfinite(tot) or not (95 <= tot <= 105):
            skip["texture"] += 1
            continue
        sand, silt, clay = 100 * frac / tot
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
        top, bottom = (float(x) for x in s["depth"].split("-"))
        intact = s["condition"].lower().startswith("intact")
        rows.append(dict(
            layer_id=layer_id(s["site"], s["depth"], s["condition"]),
            profile_id=f"BW_{s['site']}",
            texture_class=st.usda_class(sand, silt, clay),
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=sand, silt=silt, clay=clay, ksat_cmh=np.nan,
            depth_cm=(top + bottom) / 2, lat=loc.loc[s["site"], "Lat_WGS84"],
            lon=loc.loc[s["site"], "Long_WGS84"], source_db="CSIRO_Boorowa",
            oc=t["TOC%"], om=np.nan, porosity=np.nan, bd=s["bd"],
            rmse=rmse, n_points=len(h), **fm,
            sample_type="undisturbed" if intact else "disturbed",
            sample_type_source="CSIRO Boorowa: " + ("intact core" if intact
                                                    else "repacked sample")))
    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"wrote {OUT}: {len(out)} of {sum(skip.values()) + len(out)} samples; "
          f"dropped {skip}")
    print(out.groupby(["sample_type", "texture_class"]).size().to_string())


if __name__ == "__main__":
    main()
