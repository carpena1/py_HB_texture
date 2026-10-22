"""Australian intact-core retention curves from the CSIRO National Soil Site
Database (NatSoil), with USDA texture.

CSIRO (2026). CSIRO National Soil Site Database. doi:10.25919/sg5w-ny87,
CC BY 4.0. Public API: https://soils.csiro.au/data/natsoil-api/

Why: GSHP's "Australian_database" (768 layers: 653 CSIRO_2020, 115 from
Forrest et al. 1985) carries these same curves, but
  * the CSIRO_2020 layers (NatSoil projects Morphology, RAALS, CSIRO_LW,
    SSM, DD, CL, GYC) come without particle size, so prepare_gshp.py had to
    drop all but 4 of them;
  * the Forrest layers carry the Australian particle-size limits (silt
    2-20 um, sand 20-2000 um) as if they were USDA (silt 2-50 um): for 73 of
    114 the "USDA" clay and silt equal NatSoil's Coventry-Fett pipette
    values. Read that way, about 29 % of soils get the wrong USDA class;
  * the Forrest coordinates kept their decimals but lost the whole degrees
    (every site was put between 36-37 S and 142-143 E; NatSoil places the
    same sites from 25.7 to 36.7 S and 138 to 151 E).
This table rebuilds those soils from NatSoil instead, and prepare_gshp.py
leaves GSHP's Australian rows out so that nothing is counted twice.

Sites: the NatSoil sites within 100 m of a GSHP CSIRO_2020 layer, plus every
site of NatSoil's Forrest project (GSHP's Forrest coordinates cannot be
matched). download() fetches them into data/natsoil_raw/ (git-ignored), one
request per site with horizons, samples and lab results.

Curves (one per replicate core):
  * Intact cores on suction plates, NatSoil codes P3B3VL{a,b,c,e}<bar>
    (76, 73x75, 98 and 60x34 mm cores), saturation (h = 0) to 0.66 bar,
    m3/m3. Each replicate (labr_no) of a sample is a core of its own, with
    its own bulk density and, where measured, Ks -- as GSHP lists them.
  * The dry end on the same sample: pressure plates at 1-15 bar on
    disturbed (P3B2VL_), <2 mm (P3B1VL_) or small 48x15 mm intact cores
    (P3B3VLd). NatSoil records these once per sample, so their mean is
    attached to every core of the sample (GSHP attached it to the first
    replicate only, leaving the second without a dry end).
  * Potentials are read from the method description ("0.05 BAR ..."),
    1 bar = 100 kPa. Points beyond 1500 kPa are dropped, as elsewhere.
  * Kept with at least 5 points, one of them at 1 kPa or wetter, and a van
    Genuchten fit within 0.03 m3/m3. The wet-end rule is stricter than the
    10 kPa used elsewhere because the old Queensland sets (method VLe) start
    at 4 kPa and their fitted saturation then runs to 0.8-1.0 m3/m3, above
    the core's porosity; fits whose saturation exceeds 1 - BD/2.65 + 0.03
    are dropped for the same reason.

Texture (per horizon, from whichever sample carries it), USDA limits, method
recorded in texture_method:
  * sedigraph: cumulative % passing at 2, 31.2 and 53 um; the 50 um point is
    interpolated in log size between 31.2 and 53 um (a 0.5-point step).
  * pipette+sieve: Coventry-Fett pipette clay (<2 um) and silt (2-20 um)
    plus the wet-sieved 20-63 um fraction, of which the log-size share
    ln(50/20)/ln(63/20) = 0.80 is taken as 20-50 um silt.
  * pipette, interpolated (flagged): pipette clay, silt, fine sand
    (20-200 um) and coarse sand only (Forrest; and the "not recorded"
    fractions of the old Queensland surveys). The 50 um point is
    interpolated in log size between 20 and 200 um. Checked on the 188
    NatSoil horizons with a full Sedigraph curve: bias +0.6 points of silt,
    mean error 4.1, and the same USDA class as the measured split for 94 %
    of horizons (71 % if the Australian silt is read as USDA silt).
  * no particle size: the horizon is left out.

Other fields: bulk density P3A1 of the core (sample mean otherwise); Ks
P4_sat_McK (McKenzie & Jacquier 1996), P4_sat_For (Forrest et al. 1985) or
P4_sat_LOV (Loveday falling head), mm/h -> cm/h, of the core where it is
numbered so; organic carbon by Dumas combustion (6B3, 6B2) or Walkley-Black
(6A1); depth = mid-point of the sample; coordinates GDA94 (= WGS84 at this
scale), the original being AGD66. sample_type undisturbed: the wet end is on
intact cores.

Two tables (2026-10-21, project decision). All 552 cores together read
poorly as a new source (18 % exact class; always answering clay scores 41 %)
and moved a new source's exact class for the rest of the reference by
-1.6 points (p = 0.052, verify_join_effect.py), while Ks improved slightly
(-0.013 dex, p < 0.001). The misses come mostly from the Sedigraph projects
(Morphology, RAALS, CSIRO_LW: 13-29 % exact with their own sites in the
reference and depth and bulk density, against 57 % for Forrest and 30 % for
SSM), whose curves read about 10 points less clayey than measured. So only
Forrest and SSM go to natsoil_reference.csv, a candidate for the default;
the Sedigraph projects go to natsoil_other_reference.csv, a test set.
"""
import glob
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

import free_m
import swcc_texture as st

API = "https://soils.csiro.au/data/natsoil-api"
RAW = os.path.join("data", "natsoil_raw")
GSHP_RAW = os.path.join("data", "WRC_dataset_surya_et_al_2021_final.csv")
OUT = os.path.join("data", "natsoil_reference.csv")              # joins the default
OUT_OTHER = os.path.join("data", "natsoil_other_reference.csv")  # test set
DEFAULT_PROJECTS = ("Forrest", "SSM")
MIN_POINTS = 5
MAX_RMSE = 0.03
WET_KPA = 1.0              # the wettest point must be at or below this
POROSITY_SLACK = 0.03      # fitted theta_s may exceed 1 - BD/2.65 by this much
MAX_KPA = 1500.0
CORE = re.compile(r"^P3B3VL[abce]")          # intact cores on suction plates
DRY = re.compile(r"^(P3B2VL_|P3B1VL_|P3B3VLd)")   # pressure plates, once per sample
KS = ("P4_sat_McK", "P4_sat_For", "P4_sat_LOV")
OC = ("6B3", "6B2", "6A1")
LOG = np.log


def download(match_km=0.1):
    """Fetch the NatSoil sites behind GSHP's Australian layers into RAW."""
    from scipy.spatial import cKDTree
    os.makedirs(RAW, exist_ok=True)
    meta_fn = os.path.join(RAW, "site_meta.json")
    if not os.path.exists(meta_fn):
        urllib.request.urlretrieve(f"{API}/site-meta-data", meta_fn)
    rows = []
    for ag in json.load(open(meta_fn)):
        for p in ag.get("zProjects", []):
            for s in p.get("zSites", []):
                for loc in s.get("locations", []) or [{}]:
                    rows.append(dict(agency=ag["agency_code"], proj=p["proj_code"], s_id=str(s["s_id"]),
                                     lat=loc.get("latitude"), lon=loc.get("longitude")))
    ns = pd.DataFrame(rows).dropna(subset=["lat", "lon"]).drop_duplicates(["agency", "proj", "s_id"])
    g = pd.read_csv(GSHP_RAW, low_memory=False, encoding="latin-1").drop_duplicates("layer_id")
    g = g[g.reference == "CSIRO_2020"]
    sites = g.groupby(g.layer_id.str.split().str[0]).agg(
        lat=("latitude_decimal_degrees", "first"), lon=("longitude_decimal_degrees", "first"))
    k = np.cos(np.radians(-30.0))
    dist, idx = cKDTree(np.c_[ns.lat, ns.lon * k]).query(np.c_[sites.lat, sites.lon * k])
    want = pd.concat([ns.iloc[idx[dist * 111.0 < match_km]], ns[ns.proj == "Forrest"]])
    want = want[["agency", "proj", "s_id"]].drop_duplicates()
    for r in want.itertuples():
        fn = os.path.join(RAW, f"{r.agency}+{r.proj}+{r.s_id}.json".replace("/", "_"))
        if os.path.exists(fn) and os.path.getsize(fn) > 10:
            continue
        q = urllib.parse.urlencode(dict(agency_code=r.agency, proj_code=r.proj, s_id=r.s_id,
                                        horizons="true", samples="true", labResults="true"))
        open(fn, "wb").write(urllib.request.urlopen(f"{API}/sites?{q}", timeout=120).read())
        time.sleep(0.5)
    print(f"NatSoil sites in {RAW}: {len(want)}")


def kpa(desc):
    """Matric potential in kPa from a NatSoil method description."""
    if desc.lower().startswith("saturated"):
        return 0.0
    m = re.match(r"\s*([\d.]+)\s*BAR", desc)
    return float(m.group(1)) * 100.0 if m else None


def texture(v):
    """(sand, silt, clay, method) on USDA limits from a horizon's lab values."""
    g = lambda k: float(np.mean(v[k])) if k in v else np.nan
    if all(k in v for k in ("P10_S_2", "P10_S_31.2", "P10_S_53")):
        s50 = g("P10_S_31.2") + (g("P10_S_53") - g("P10_S_31.2")) * LOG(50 / 31.2) / LOG(53 / 31.2)
        clay = g("P10_S_2")
        return 100.0 - s50, s50 - clay, clay, "sedigraph"
    if all(k in v for k in ("P10_CF_C", "P10_CF_Z", "PWS20-63")):
        clay = g("P10_CF_C")
        silt = g("P10_CF_Z") + LOG(50 / 20) / LOG(63 / 20) * g("PWS20-63")
        return 100.0 - clay - silt, silt, clay, "pipette+sieve"
    for pre in ("P10_CF_", "P10_NR_"):
        if all(pre + k in v for k in ("C", "Z", "FS", "CS")):
            c, z, fs, cs = (g(pre + k) for k in ("C", "Z", "FS", "CS"))
            tot = c + z + fs + cs
            if not 90.0 <= tot <= 110.0:
                return (np.nan,) * 3 + (None,)
            c, z, fs, cs = (100.0 * x / tot for x in (c, z, fs, cs))
            silt = z + LOG(50 / 20) / LOG(200 / 20) * fs
            return 100.0 - c - silt, silt, c, "pipette, interpolated 20-200 um"
    return (np.nan,) * 3 + (None,)


def cores():
    """One record per replicate core: points, texture and covariates."""
    out = []
    for fn in sorted(glob.glob(os.path.join(RAW, "*+*.json"))):
        for site in json.load(open(fn)):
            for o in site.get("z_observations", []):
                lat = o.get("o_latitude_GDA94") or o.get("o_latitude")
                lon = o.get("o_longitude_GDA94") or o.get("o_longitude")
                for h in o.get("z_horizons", []):
                    hv = {}
                    for s in h.get("z_samples", []) or []:
                        for r in s.get("z_lab_results") or []:
                            if r["labr_value"] is not None:
                                hv.setdefault(r["labm_code"], []).append(float(r["labr_value"]))
                    sand, silt, clay, tmeth = texture(hv)
                    oc = next((float(np.mean(hv[k])) for k in OC if k in hv), np.nan)
                    for s in h.get("z_samples", []) or []:
                        lr = [r for r in (s.get("z_lab_results") or []) if r["labr_value"] is not None]
                        dry = {}
                        for r in lr:
                            if DRY.match(r["labm_code"]) and r.get("labm_units") == "m3/m3":
                                p = kpa(r["labm_desc"] or "")
                                if p is not None and p <= MAX_KPA:
                                    dry.setdefault(p, []).append(float(r["labr_value"]))
                        reps = sorted({r["labr_no"] for r in lr if CORE.match(r["labm_code"])})
                        for rep in reps:
                            pts = {}
                            for r in lr:
                                if r["labr_no"] == rep and CORE.match(r["labm_code"]) \
                                        and r.get("labm_units") == "m3/m3":
                                    p = kpa(r["labm_desc"] or "")
                                    if p is not None:
                                        pts.setdefault(p, []).append(float(r["labr_value"]))
                            for p, vals in dry.items():
                                pts.setdefault(p, vals)
                            pick = lambda codes: next(
                                (float(r["labr_value"]) for c in codes for r in lr
                                 if r["labm_code"] == c and r["labr_no"] == rep), None)
                            bd = pick(("P3A1",))
                            if bd is None:
                                bd = float(np.mean([float(r["labr_value"]) for r in lr if r["labm_code"] == "P3A1"])) \
                                    if any(r["labm_code"] == "P3A1" for r in lr) else np.nan
                            ks = pick(KS)
                            top, bot = s.get("samp_upper_depth"), s.get("samp_lower_depth")
                            depth = (100.0 * (float(top) + float(bot)) / 2.0
                                     if top is not None and bot is not None else np.nan)
                            out.append(dict(
                                layer_id=f"NS_{site['proj_code']}_{site['s_id']}_{o['o_id']}_{h['h_no']}_{s['samp_no']}_{rep}",
                                profile_id=f"NS_{site['agency_code']}_{site['proj_code']}_{site['s_id']}",
                                project=site["proj_code"], points=pts, sand=sand, silt=silt, clay=clay,
                                texture_method=tmeth, oc=oc, bd=bd,
                                ksat_cmh=ks / 10.0 if ks is not None else np.nan,
                                depth_cm=depth, lat=float(lat) if lat else np.nan,
                                lon=float(lon) if lon else np.nan))
    return out


def measured_points():
    """{layer_id: (h kPa, theta)} for every kept core, for tests that refit."""
    out = {}
    for c in cores():
        h = np.array(sorted(c["points"]))
        th = np.array([np.mean(c["points"][p]) for p in h])
        out[c["layer_id"]] = (h, th)
    return out


def main():
    if "--download" in sys.argv or not glob.glob(os.path.join(RAW, "*+*.json")):
        download()
    rows, skip = [], dict(points=0, texture=0, fit=0, rmse=0, porosity=0)
    for c in cores():
        h = np.array(sorted(c["points"]))
        th = np.array([np.mean(c["points"][p]) for p in h])
        if len(h) < MIN_POINTS or not (h <= WET_KPA).any():
            skip["points"] += 1
            continue
        if c["texture_method"] is None:
            skip["texture"] += 1
            continue
        try:
            popt, _ = st.fit_vg(h, th)
        except Exception:
            skip["fit"] += 1
            continue
        tr, ts, la, ln1 = popt
        alpha, n = 10.0 ** la, 1.0 + 10.0 ** ln1
        rmse = float(np.sqrt(np.mean((st.vg_theta(h, tr, ts, alpha, n) - th) ** 2)))
        if not np.isfinite(rmse) or rmse > MAX_RMSE or n <= 1.0:
            skip["rmse"] += 1
            continue
        if np.isfinite(c["bd"]) and ts > 1.0 - c["bd"] / 2.65 + POROSITY_SLACK:
            skip["porosity"] += 1
            continue
        fm = free_m.fit_free_m(h, th, fallback=free_m.mualem_fallback(tr, ts, alpha, n))
        sand, silt, clay = c["sand"], c["silt"], c["clay"]
        rows.append(dict(
            layer_id=c["layer_id"], profile_id=c["profile_id"],
            texture_class=st.usda_class(sand, silt, clay),
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=round(sand, 2), silt=round(silt, 2), clay=round(clay, 2),
            ksat_cmh=c["ksat_cmh"], depth_cm=c["depth_cm"], lat=c["lat"], lon=c["lon"],
            source_db="NatSoil", oc=c["oc"], om=1.724 * c["oc"] if np.isfinite(c["oc"]) else np.nan,
            bd=c["bd"], rmse=rmse, n_points=len(h), **fm,
            sample_type="undisturbed",
            sample_type_source="NatSoil: intact cores on suction plates; 1-15 bar "
                               "pressure plates on the same sample",
            texture_method=c["texture_method"], natsoil_project=c["project"]))
    out = pd.DataFrame(rows)
    keep = out.natsoil_project.isin(DEFAULT_PROJECTS)
    out[keep].to_csv(OUT, index=False)
    out[~keep].to_csv(OUT_OTHER, index=False)
    print(f"{OUT}: {keep.sum()} cores ({', '.join(DEFAULT_PROJECTS)}); {OUT_OTHER}: {(~keep).sum()}")
    print(f"cores: {len(out)} kept   (skipped: {skip})")
    print(f"  sites {out.profile_id.nunique()}; projects {out.natsoil_project.value_counts().to_dict()}")
    print(f"  texture method {out.texture_method.value_counts().to_dict()}")
    print(f"  with Ks {out.ksat_cmh.notna().sum()}, bd {out.bd.notna().sum()}, oc {out.oc.notna().sum()}, "
          f"coordinates {out.lat.notna().sum()}")
    print(f"  median fit RMSE {out.rmse.median():.4f}; points per curve {out.n_points.median():.0f}")
    print(out.texture_class.value_counts().to_string())


if __name__ == "__main__":
    main()
