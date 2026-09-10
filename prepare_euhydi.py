"""One-time preprocessing of EU-HYDI, the European Hydropedological Data
Inventory (Weynants et al. 2013, EUR 26053 EN, doi:10.2788/5936).

Reads data/euhydi_raw/HYDI-v1.accdb (needs mdbtools) and writes
data/euhydi_reference.csv.

RESTRICTED DATA. EU-HYDI is available to consortium members only, and the
intellectual property stays with the individual data providers. Both the raw
database and this derived table are gitignored; neither may be committed to
the public repository until redistribution terms are agreed. The table is
part of the default `merged` reference wherever it has been built; a checkout
without it falls back to the public tables (see
swcc_texture.RESTRICTED_TABLES).

EU-HYDI includes the sample-level HYPRES data (Wosten_HYPRES,
Schindler_HYPRES, Hennings_HYPRES, ...) that the ESDB distribution never
released. Each contributor is kept as its own source_db, so source-blocked
validation treats them as separate laboratories.

Choices, each matching an existing preparation script:
  * Retention points flagged FLAG = FALSE ("unreliable") are dropped, and a
    sample is dropped if any water content exceeds MAX_THETA.
  * Points are thinned to one median per PF_BIN of log suction, as in
    prepare_hohenbrink.py: some contributors ship evaporation-method curves
    with up to ~1,200 points, which would otherwise outweigh the dry end.
  * Curves whose wettest point is drier than TRUNCATED_KPA get the same
    saturation anchor as prepare_kssl.py, theta(0) = 0.95 x porosity. This
    affects about 3 % of samples; 97 % already reach 6 kPa or wetter.
  * Texture uses the USDA-system fractions in PSD_EST, accepting "measured"
    and "interpolated". Most European laboratories use a 63 um silt/sand
    break, and EU-HYDI harmonised to USDA's 50 um by interpolating the
    measured grain-size curve. Where PSD_EST holds no values -- all of the
    Polish data, coded "PSD measured up to 1000 microns" but left empty --
    the raw size classes in PSIZE are summed instead, and only when they
    break at exactly 2 and 50 um, so no interpolation of our own is involved.
    Sand then stops at 1 mm rather than 2 mm, which renormalising absorbs.
  * Not usable, by construction rather than by filter: Bilas (Greece) and
    Anaya (Spain) report a median of 2 retention points per sample, Hennings
    exactly 4, and Patyka (Ukraine) none -- too few to fit 4 parameters.
  * Ksat is K at h = 0 from saturated-conductivity methods only (codes
    800-819 laboratory cores, 835 column, 850-869 in situ). Evaporation,
    crust and hot-air methods reach h = 0 only by extrapolation.

Units:
    HEAD   cm suction      ->  h_kPa = HEAD * 0.0980665
    THETA  cm3/cm3
    COND   cm/day          ->  ksat_cmh = COND / 24
    POR    volume %        ->  porosity = POR / 100
    OC     mass %
"""

import io
import os
import subprocess

import numpy as np
import pandas as pd

import free_m
import swcc_texture as st

RAW = os.path.join("data", "euhydi_raw", "HYDI-v1.accdb")
OUT = os.path.join("data", "euhydi_reference.csv")
CM_TO_KPA = 0.0980665
PF_BIN = 0.25
MIN_POINTS = 5
MAX_RMSE = 0.03
MAX_THETA = 0.95
TRUNCATED_KPA = 6.0
SATURATION_FRACTION = 0.95
PARTICLE_DENSITY = 2.65
PSD_OK = {"measured", "interpolated"}
SENTINELS = [-999, -999.0, -9999, -9999.0]


def is_sat_method(code):
    return (800 <= code <= 819) or code == 835 or (850 <= code <= 869)


def table(name):
    """One Access table as a DataFrame, with -999 sentinels as NaN."""
    txt = subprocess.run(["mdb-export", RAW, name], capture_output=True,
                         text=True, check=True).stdout
    df = pd.read_csv(io.StringIO(txt), low_memory=False)
    return df.replace(SENTINELS, np.nan)


def thin(h_cm, theta):
    """One median point per PF_BIN of log suction; h = 0 points kept as one."""
    sat = h_cm == 0
    out_h, out_t = [], []
    if sat.any():
        out_h.append(0.0)
        out_t.append(float(np.median(theta[sat])))
    if (~sat).any():
        g = pd.DataFrame({"pf": np.log10(h_cm[~sat]), "th": theta[~sat]})
        g = g.groupby(np.floor(g.pf / PF_BIN)).median()
        out_h += list(10.0 ** g.pf.to_numpy())
        out_t += list(g.th.to_numpy())
    return np.array(out_h), np.array(out_t)


def usda_from_size_classes(psize):
    """USDA (sand, silt, clay) % from raw PSIZE classes, for samples whose
    classes break at exactly 2 and 50 um. P_PERCENT is the mass % between a
    class's P_SIZE and the next smaller P_SIZE of the same sample."""
    out = {}
    for sid, g in psize.dropna(subset=["P_SIZE", "P_PERCENT"]).groupby("SAMPLE_ID"):
        if not {2.0, 50.0} <= set(g.P_SIZE):
            continue
        out[sid] = (g.P_PERCENT[(g.P_SIZE > 50) & (g.P_SIZE <= 2000)].sum(),
                    g.P_PERCENT[(g.P_SIZE > 2) & (g.P_SIZE <= 50)].sum(),
                    g.P_PERCENT[g.P_SIZE <= 2].sum())
    return out


def main():
    if not os.path.exists(RAW):
        raise SystemExit(f"{RAW} not found -- EU-HYDI is consortium-restricted "
                         f"and is not distributed with this repository.")
    ret = table("RET")
    basic = table("BASIC").set_index("SAMPLE_ID")
    psd = table("PSD_EST").set_index("SAMPLE_ID")
    gen = table("GENERAL").set_index("PROFILE_ID")
    chem = table("CHEMICAL").set_index("SAMPLE_ID")
    cond = table("COND")
    psize = table("PSIZE")
    for c in ("P_SIZE", "P_PERCENT"):
        psize[c] = pd.to_numeric(psize[c], errors="coerce")
    classes = usda_from_size_classes(psize)

    for c in ("HEAD", "THETA"):
        ret[c] = pd.to_numeric(ret[c], errors="coerce")
    ret = ret[ret.FLAG.astype(str) != "False"]
    ret = ret[ret.HEAD.notna() & ret.THETA.notna() & (ret.HEAD >= 0)]

    for c in ("IND_VALUE", "VALUE", "COND", "COND_M"):
        cond[c] = pd.to_numeric(cond[c], errors="coerce")
    ks = cond[(cond.IND_VALUE == 1) & (cond.VALUE == 0) & (cond.COND > 0)
              & cond.COND_M.fillna(-1).map(is_sat_method)]
    ksat = ks.groupby("SAMPLE_ID").COND.median() / 24.0

    rows = []
    skip = dict(points=0, theta=0, texture=0, porosity=0, fit=0, rmse=0)
    n_anchor = 0
    for sid, g in ret.groupby("SAMPLE_ID"):
        h_cm = g.HEAD.to_numpy(float)
        th = g.THETA.to_numpy(float)
        if th.max() > MAX_THETA or th.min() <= 0:
            skip["theta"] += 1
            continue
        h_cm, th = thin(h_cm, th)
        if len(h_cm) < MIN_POINTS:
            skip["points"] += 1
            continue
        h = h_cm * CM_TO_KPA

        frac, code = None, None
        if sid in psd.index:
            p = psd.loc[sid]
            if isinstance(p, pd.DataFrame):
                p = p.iloc[0]
            if {p.USCLAY_CODE, p.USSILT_CODE, p.USSAND_CODE} <= PSD_OK:
                frac = (float(p.USSAND), float(p.USSILT), float(p.USCLAY))
                code = p.USSILT_CODE
        if (frac is None or not np.isfinite(sum(frac))) and sid in classes:
            frac, code = classes[sid], "size classes"
        tot = sum(frac) if frac is not None else np.nan
        if not np.isfinite(tot) or not (95 <= tot <= 105):
            skip["texture"] += 1
            continue
        sand, silt, clay = (100 * x / tot for x in frac)

        b = basic.loc[sid] if sid in basic.index else None
        if isinstance(b, pd.DataFrame):
            b = b.iloc[0]
        bd = float(b.BD) if b is not None and pd.notna(b.BD) else np.nan
        por = (float(b.POR) / 100.0 if b is not None and pd.notna(b.POR)
               else 1.0 - bd / PARTICLE_DENSITY)

        n_measured = len(h)
        if h.min() > TRUNCATED_KPA:
            # Same wet-end fix as prepare_kssl.py: nothing constrains
            # saturation, so anchor it, floored just above the wettest point.
            if not np.isfinite(por) or not (0.2 < por < 0.9):
                skip["porosity"] += 1
                continue
            anchor = max(por * SATURATION_FRACTION, th.max() * 1.05)
            if anchor > MAX_THETA:
                skip["porosity"] += 1
                continue
            h = np.concatenate([[0.0], h])
            th = np.concatenate([[anchor], th])
            n_anchor += 1

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
        fm = free_m.fit_free_m(h, th,
                               fallback=free_m.mualem_fallback(tr, ts, alpha, n))

        pid = g.PROFILE_ID.iloc[0]
        gi = gen.loc[pid] if pid in gen.index else None
        lat = float(gi.Y_WGS84) if gi is not None and pd.notna(gi.Y_WGS84) else np.nan
        lon = float(gi.X_WGS84) if gi is not None and pd.notna(gi.X_WGS84) else np.nan
        src = gi.SOURCE if gi is not None and pd.notna(gi.SOURCE) else g.SOURCE.iloc[0]
        top = float(b.SAMPLE_DEP_TOP) if b is not None and pd.notna(b.SAMPLE_DEP_TOP) else np.nan
        bot = float(b.SAMPLE_DEP_BOT) if b is not None and pd.notna(b.SAMPLE_DEP_BOT) else np.nan
        oc = np.nan
        if sid in chem.index:
            c = chem.loc[sid]
            c = c.iloc[0] if isinstance(c, pd.DataFrame) else c
            oc = float(c.OC) if pd.notna(c.OC) and 0 <= float(c.OC) < 60 else np.nan

        rows.append(dict(
            layer_id=f"EUHYDI_{sid}", profile_id=f"EUHYDI_{pid}",
            texture_class=st.usda_class(sand, silt, clay),
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=sand, silt=silt, clay=clay,
            ksat_cmh=ksat.get(sid, np.nan),
            depth_cm=(top + bot) / 2.0 if np.isfinite(top + bot) else top,
            lat=lat, lon=lon, source_db=f"EUHYDI_{src}",
            oc=oc, porosity=por, rmse=rmse, n_points=n_measured,
            country=gi.ISO_COUNTRY if gi is not None else np.nan,
            psd_code=code, **fm))

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"samples with retention data: {ret.SAMPLE_ID.nunique()}")
    print(f"reference layers: {len(out)}   (skipped: {skip})")
    print(f"  wet-end anchored:  {n_anchor}")
    print(f"  with ksat:         {out.ksat_cmh.notna().sum()}   "
          f"(median {out.ksat_cmh.median():.2f} cm/h)")
    print(f"  with coordinates:  {out.lat.notna().sum()}")
    print(f"  with oc:           {out.oc.notna().sum()}")
    print(f"  with depth:        {out.depth_cm.notna().sum()}")
    print(f"  profiles:          {out.profile_id.nunique()}   "
          f"contributors: {out.source_db.nunique()}")
    print(f"  texture: {out.psd_code.value_counts().to_dict()}")
    print(f"  median fit RMSE:   {out.rmse.median():.4f}")

    # Unit check: class-median alpha and n against GSHP, as prepare_kssl.py.
    gs = pd.read_csv(os.path.join("data", "gshp_reference.csv"))
    print("\nunit check -- class-median alpha [kPa^-1] / n, EU-HYDI vs GSHP:")
    for cl in ("sandy loam", "loam", "silt loam", "clay loam", "clay"):
        a, b = out[out.texture_class == cl], gs[gs.texture_class == cl]
        print(f"  {cl:<12s} alpha {a.alpha_kpa.median():7.3f} vs "
              f"{b.alpha_kpa.median():7.3f}   n {a.n.median():5.2f} vs "
              f"{b.n.median():5.2f}   ({len(a)} vs {len(b)} soils)")
    print()
    print(out.texture_class.value_counts().to_string())
    print()
    print(out.country.value_counts().to_string())


if __name__ == "__main__":
    main()
