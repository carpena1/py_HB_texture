"""One-time preprocessing of the New Jersey soils of Soil Survey Investigations
Report No. 26.

USDA Soil Conservation Service (1974). Soil survey laboratory data and
descriptions for some soils of New Jersey. Soil Survey Investigations
Report No. 26, in cooperation with the New Jersey Agricultural Experiment
Station, Rutgers University. A US government publication (HathiTrust scan
uc1.d0005445648).

The same soils, and the reason they are here: Arya, L.M., Richter, J.C. and
Davidson, S.A. (1982). A comparison of soil moisture characteristics
predicted by the Arya-Paris model with laboratory-measured data. AgRISTARS
report SM-L2-04247 (JSC-17820), NASA Johnson Space Center -- the Arya-Paris
model's own large test set, 181 of these horizons.

Reads, all under data/arya_nj_raw/ (git-ignored):
  * ssir26/batch_*.jsonl -- the report's data sheets transcribed from the page
    images (the scan's text layer is unusable), one pedon per line.
  * transcribed/batch_*.jsonl -- the NASA report's appendix B, transcribed the
    same way, used here only as an independent check of the curves: its
    volumetric water contents are Rutgers' weight percentages times bulk
    density, copied and printed a second time.
and the NCSS lab database (data/kssl_raw/NCSSLabDataMartSQLite.sqlite3) for a
check of the texture, which is the same Beltsville pipette analysis. Writes
data/nj_ssir26_reference.csv.

  * 46 pedons with Rutgers analyses, 1956-61, New Jersey Coastal Plain.
  * Retention (Rutgers): 0.02, 0.06, 0.1 and 0.33 bar on a tension table and
    1 bar in a pressure membrane on saturated 3-inch cores; 2, 6 and 15 bar on
    crushed samples. Weight percent, converted with the core's bulk density.
    An undisturbed wet end with a disturbed dry end, counted as undisturbed
    as elsewhere in the project. The splice shows: in some horizons 2 bar
    holds more water than 1 bar. Points are kept as measured.
  * Saturation: the total porosity, "water retained at 0 bars suction:
    saturated", measured on the same cores. It is added at h = 0 only when
    the sheet's non-capillary + capillary = total holds (within 0.5 point);
    where a typo breaks it, the point is left out and the curve fitted
    without it.
  * Ks: constant head (1 inch of water) on the saturated cores, in/hr ->
    cm/h: an intact-core Ks, which Zanjanrood's repacked samples are not.
  * Texture: the SCS pipette sand, silt and clay of the same sheet.
    Rows whose three do not sum to 95-105 are left out (a few printed
    typos), the rest renormalised.
"""

import glob
import json
import os
import re
import sqlite3

import numpy as np
import pandas as pd

import free_m
import swcc_texture as st

RAW = os.path.join("data", "arya_nj_raw")
NCSS = os.path.join("data", "kssl_raw", "NCSSLabDataMartSQLite.sqlite3")
OUT = os.path.join("data", "nj_ssir26_reference.csv")
BARS = ["0.02", "0.06", "0.1", "0.33", "1", "2", "6", "15"]
KPA_PER_BAR = 100.0
CM_PER_IN = 2.54
OM_PER_OC = 1.724
MIN_POINTS = 5
MAX_RMSE = 0.03
POROSITY_TOL = 0.5       # points: non-capillary + capillary against total
AW_TOL = 0.02            # in/in: available water against (FC - 15 bar) x BD

# Printed typos, found where the sheet and the NASA report disagree and
# settled by the sheet's own available water, AW = (FC - w15) x BD / 100.
CORRECTIONS = {
    ("S58NJ-5-1", "28-32", "15"): (3.1, "sheet 7.1 makes the curve rise; NASA "
                                   "1982 has 3.1 (0.051 / 1.64), and AW 0.115 "
                                   "needs 3.1-3.5"),
    ("S57NJ-12-3", "38-50", "15"): (3.6, "sheet 4.7; NASA 1982 has 3.6 (0.067 / "
                                    "1.85), and AW 0.14 needs 3.6"),
}


def num(v):
    return float(v) if isinstance(v, (int, float)) else np.nan


def depth_cm(s):
    """'0-10' inches -> (top, bottom) cm; an open '30' -> (76.2, nan)."""
    parts = [p for p in re.split(r"\s*-\s*", str(s)) if p]
    top = float(parts[0]) * CM_PER_IN
    bot = float(parts[1]) * CM_PER_IN if len(parts) > 1 else np.nan
    return top, bot


def read_sheets():
    pages = []
    for path in sorted(glob.glob(os.path.join(RAW, "ssir26", "batch_*.jsonl"))):
        pages += [json.loads(line) for line in open(path) if line.strip()]
    pages = [p for p in pages if p["horizons"]]
    for p in pages:
        for h in p["horizons"]:
            for b in BARS:
                fix = CORRECTIONS.get((p["soil_no"], str(h["depth_in"]), b))
                if fix:
                    h["w_bar"][b] = fix[0]
    return pages


def aw_check(sheets):
    """Rows whose printed available water disagrees with their own field
    capacity, 15-bar water and bulk density."""
    bad = []
    for p in sheets:
        for h in p["horizons"]:
            aw, fc, bd = num(h.get("available_water")), num(h.get("field_capacity")), num(h.get("bd"))
            w15 = num((h.get("w_bar") or {}).get("15"))
            if np.isfinite(aw) and np.isfinite(fc) and np.isfinite(bd) and np.isfinite(w15):
                calc = (fc - w15) * bd / 100.0
                if abs(calc - aw) > AW_TOL:
                    bad.append((p["soil_no"], h["depth_in"], aw, round(calc, 3)))
    return bad


def layer_id(soil_no, depth_in):
    return f"NJ26_{soil_no}_{str(depth_in).replace(' ', '')}"


def points(h):
    """(h kPa, theta) from one sheet row, and whether saturation was added."""
    bd = num(h.get("bd"))
    w = h.get("w_bar") or {}
    hk = [float(b) * KPA_PER_BAR for b in BARS if np.isfinite(num(w.get(b)))]
    th = [num(w[b]) * bd / 100.0 for b in BARS if np.isfinite(num(w.get(b)))]
    total, nc, cp = (num(h.get(k)) for k in ("total_porosity", "noncap", "cap"))
    sat = (np.isfinite(total) and np.isfinite(nc) and np.isfinite(cp)
           and abs(nc + cp - total) <= POROSITY_TOL)
    if sat:
        hk, th = [0.0] + hk, [total / 100.0] + th
    return np.array(hk), np.array(th), sat


def measured_points():
    """{layer_id: (h kPa, theta)} for every usable row, for tests that refit."""
    out = {}
    for p in read_sheets():
        for h in p["horizons"]:
            if np.isfinite(num(h.get("bd"))):
                hk, th, _ = points(h)
                out[layer_id(p["soil_no"], h["depth_in"])] = (hk, th)
    return out


def nasa_check(sheets):
    """Compare the NASA report's volumetric points with weight % x bulk
    density, horizon by horizon (same pedon, same depth)."""
    t4 = pd.read_csv(os.path.join(RAW, "table4_horizons.csv"), dtype=str)
    nasa = {}
    for path in glob.glob(os.path.join(RAW, "transcribed", "batch_*.jsonl")):
        for line in open(path):
            if line.strip():
                r = json.loads(line)
                nasa[int(r["horizon"])] = r
    rows = {(p["soil_no"], str(h["depth_in"]).replace(" ", "")): h
            for p in sheets for h in p["horizons"]}
    diffs, bad, n_match = [], [], 0
    for _, m in t4.iterrows():
        k = int(m.horizon)
        if k not in nasa:
            continue
        key = (m.survey_code, f"{m.top_in}-{m.bottom_in}".rstrip("-"))
        h = rows.get(key)
        if h is None:
            continue
        n_match += 1
        bd = num(h.get("bd"))
        for th_n, h_cm in nasa[k]["measured"]:
            if th_n is None:
                continue
            b = min(BARS, key=lambda x: abs(np.log(float(x) * 1019.7 / h_cm)))
            w = num((h.get("w_bar") or {}).get(b))
            if np.isfinite(w) and np.isfinite(bd):
                d = th_n - w * bd / 100.0
                diffs.append(d)
                if abs(d) > 0.011:
                    bad.append((k, m.survey_code, key[1], b, th_n, round(w * bd / 100, 3)))
        if abs(num(nasa[k]["bulk_density"]) - bd) > 0.005:
            bad.append((k, m.survey_code, key[1], "bd", nasa[k]["bulk_density"], bd))
    return n_match, np.array(diffs), bad


def ncss_check(out):
    """SCS pipette clay on the sheets against the NCSS database's, same pedon
    and overlapping depth."""
    db = sqlite3.connect(f"file:{NCSS}?mode=ro", uri=True)
    d = []
    for _, r in out.iterrows():
        yy, cty, nn = re.fullmatch(r"S(\d\d)NJ-(\d+)-(\d+)", r.soil_no).groups()
        uid = f"S19{yy}NJ{2 * int(cty) - 1:03d}{int(nn):03d}"
        q = db.execute("""SELECT l.hzn_top, l.hzn_bot, p.clay_total FROM
            lab_combine_nasis_ncss c JOIN lab_layer l USING(pedon_key)
            JOIN lab_physical_properties_vw p USING(layer_key)
            WHERE c.upedonid = ? AND p.clay_total IS NOT NULL""", (uid,)).fetchall()
        top = r.depth_top_cm
        best = min(q, key=lambda x: abs(x[0] - top), default=None)
        if best is not None and abs(best[0] - top) <= 2.0:
            d.append(r.clay - best[2])
    return np.array(d)


def main():
    sheets = read_sheets()
    rows = []
    skip = dict(texture=0, bulk_density=0, points=0, fit=0, rmse=0)
    no_sat, non_mono = 0, 0
    for p in sheets:
        for h in p["horizons"]:
            sand, silt, clay = (num(h.get(c)) for c in ("sand", "silt", "clay"))
            tot = sand + silt + clay
            if not np.isfinite(tot) or not 95 <= tot <= 105:
                skip["texture"] += 1
                continue
            sand, silt, clay = 100 * sand / tot, 100 * silt / tot, 100 * clay / tot
            bd = num(h.get("bd"))
            if not 0.5 <= bd <= 2.2:
                skip["bulk_density"] += 1
                continue
            hk, th, sat = points(h)
            if len(hk) < MIN_POINTS:
                skip["points"] += 1
                continue
            no_sat += not sat
            non_mono += bool(np.any(np.diff(th) > 0.005))
            try:
                popt, _ = st.fit_vg(hk, th)
            except Exception:
                skip["fit"] += 1
                continue
            thr, ths, la, ln1 = popt
            alpha, n = 10.0 ** la, 1.0 + 10.0 ** ln1
            rmse = float(np.sqrt(np.mean((st.vg_theta(hk, thr, ths, alpha, n) - th) ** 2)))
            if not np.isfinite(rmse) or rmse > MAX_RMSE or n <= 1.0:
                skip["rmse"] += 1
                continue
            fm = free_m.fit_free_m(hk, th,
                                   fallback=free_m.mualem_fallback(thr, ths, alpha, n))
            top, bot = depth_cm(h["depth_in"])
            oc = num(h.get("oc"))
            ks = num(h.get("ksat_in_hr")) * CM_PER_IN
            rows.append(dict(
                layer_id=layer_id(p["soil_no"], h["depth_in"]),
                profile_id=f"NJ26_{p['soil_no']}",
                texture_class=st.usda_class(sand, silt, clay),
                alpha_kpa=alpha, n=n, thetar=thr, thetas=ths,
                sand=sand, silt=silt, clay=clay,
                ksat_cmh=ks if ks > 0 else np.nan,
                depth_cm=(top + bot) / 2 if np.isfinite(bot) else top,
                lat=np.nan, lon=np.nan, source_db="SSIR26_NewJersey",
                oc=oc, om=OM_PER_OC * oc if np.isfinite(oc) else np.nan,
                porosity=num(h.get("total_porosity")) / 100.0, bd=bd,
                rmse=rmse, n_points=len(hk), **fm,
                sample_type="undisturbed",
                sample_type_source="SSIR 26 (Rutgers): 0-1 bar and Ks on "
                                   "saturated 3-inch cores, 2-15 bar crushed",
                soil_no=p["soil_no"], soil=p["soil"], county=p["county"],
                horizon=h.get("horizon"), depth_top_cm=top,
                saturation_point=sat, pdf_page=p["pdf_page"]))

    out = pd.DataFrame(rows)
    out.drop(columns=["soil_no", "depth_top_cm"]).to_csv(OUT, index=False)
    print(f"sheets with Rutgers data: {len(sheets)}; horizon rows "
          f"{sum(len(p['horizons']) for p in sheets)}")
    print(f"reference layers: {len(out)}   (skipped: {skip})")
    print(f"  pedons {out.profile_id.nunique()}; with Ks {out.ksat_cmh.notna().sum()} "
          f"(median {out.ksat_cmh.median():.2f} cm/h); saturation point left out "
          f"for {no_sat}; curves with a rise across the core/crushed splice "
          f"{non_mono}")
    print(f"  median fit RMSE {out.rmse.median():.4f}; bulk density median "
          f"{out.bd.median():.2f}")
    n_match, diffs, bad = nasa_check(sheets)
    if len(diffs):
        print(f"  NASA report check: {n_match} horizons matched, {len(diffs)} points; "
              f"|theta_NASA - w x BD| median {np.median(np.abs(diffs)):.4f}, "
              f"{np.mean(np.abs(diffs) <= 0.011) * 100:.1f} % within 0.011")
        for b in bad:
            print(f"    mismatch: NASA horizon {b[0]} {b[1]} {b[2]} in, {b[3]}: "
                  f"NASA {b[4]} vs sheet {b[5]}")
    aw = aw_check(sheets)
    print(f"  available-water check: {len(aw)} rows off by more than {AW_TOL} in/in "
          f"(after {len(CORRECTIONS)} corrections)")
    for b in aw:
        print(f"    {b[0]} {b[1]} in: printed AW {b[2]}, from FC, 15 bar and BD {b[3]}")
    dc = ncss_check(out)
    if len(dc):
        print(f"  NCSS check: clay against the database, {len(dc)} layers, "
              f"|diff| median {np.median(np.abs(dc)):.2f}, max {np.max(np.abs(dc)):.1f}")
    print(out.texture_class.value_counts().to_string())


if __name__ == "__main__":
    main()
