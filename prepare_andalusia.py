"""One-time preprocessing of K. Vanderlinden's Andalusian retention curves.

Reads data/andalusia_raw/ (shared with this project by K. Vanderlinden,
IFAPA; the site notes in each subfolder describe the files) and writes one
table per site:

  data/setenil_reference.csv   Setenil de las Bodegas (Cadiz), olive grove
  data/tomejil_reference.csv   Finca Tomejil, Carmona (Seville), vertisol
  data/donana_reference.csv    El Abalario and El Acebuche (Donana), sands
  data/lebrija_reference.csv   Sector B-XII, Lebrija (Seville), reclaimed
                               marsh clay

Answers to the project's questions are from K. Vanderlinden (email,
2026-10-02).

Setenil (Setenil/SWRC_Setenil_Rafa.xlsx):
  * 44 undisturbed rings, 11 locations x (under the canopy, inter-row) x
    (5, 15 cm). profile_id is the location; depth_cm as given.
  * Gravimetric water content at 0-500 cm on sand-kaolin suction tables
    (sheet "Final", identical to CR_columna) plus the authors' selected WP4
    dew-point readings, 3.5e3-2.2e6 cm. Volumetric = gravimetric x bulk
    density.
  * Bulk density, texture and stone content from "Datos Agrario_anillos"
    (the lab's analysis of the soil in each ring). CR_columna carries a
    density 1.0871 times larger for every ring; with it the wettest point
    would be 1.25 times the porosity in every ring (impossible), with the
    lab's value 1.01 times (median), so the lab's value is used. Confirmed:
    "Datos Agrario" is the density of the 100 cm3 rings the curves were
    measured on; the factor most likely comes from a comparison with larger
    (250 cm3) rings at some points.
  * Texture by the Laboratorio Agrario, most likely hydrometer, USDA limits
    (sand > 50 um) and carbonates removed (K. Vanderlinden).
  * The WP4 readings in "Final" were selected from WP4_analisis by visual
    inspection (the instrument is very sensitive to laboratory temperature
    and humidity); there is no stated criterion.
  * Ks: Eijkelkamp laboratory permeameter on the same rings, cm/h.
  * Organic matter from the disturbed bag next to each ring; oc = om/1.724.
  * Coordinates: UTM ED50 zone 30 (Puntos muestreo suelo_X_C.xls), one per
    location and position.

Tomejil (Tomejil/SWRCdata Tomejil.xls; Vanderlinden et al. 2021, EJSS):
  * 54 undisturbed rings (0-5 cm), 27 direct drill (DD) and 27
    conventional tillage (CT), sampled 6 October 2006. One ring per
    location, so profile_id = location. No Ks.
  * Gravimetric water content at 1-500 cm (sand and sand-kaolin boxes),
    1,000 and 3,000 cm (pressure plate) and WP4-TE readings to ~3e6 cm.
  * No per-ring bulk density. The soil is rich in smectite and its density
    runs from about 1.0 (wet) to 1.6 g/cm3 (dry), which is why the authors
    work in gravimetric water content (K. Vanderlinden). Water is made
    volumetric with a density that follows the water content, by normal
    shrinkage (the soil loses volume as it loses water): specific volume
    1/2.65 + w, density its inverse, capped at 1.75 g/cm3 (the shrinkage
    limit, w ~ 0.19); theta = w * density. Saturation then sits at density
    ~1.16 and the dry end at 1.75. This is a model, not a measured shrinkage
    curve. Tested as a new source (2026-10-02, K. Vanderlinden's question):
    a single treatment density (1.085 CT / 1.125 DD, Ordonez Fernandez et al.
    2007) reads 87 % of rings as clay with 40 % clay predicted; the
    shrinkage density reads all of them as clay with 53.6 % (measured 55.6 %).
    The cap (2026-10-12) is the median oven-dry clod density of KSSL
    Vertisols with 45-65 % clay and a 33 kPa clod density below 1.25 g/cm3
    (828 horizons, 324 pedons; 10-90 %: 1.54-1.90), the density at which
    such clods stop shrinking. It replaces 1.6, the dry end the authors
    report for Tomejil, which gives the same texture (100 % clay, 53.6 %)
    and lower dry-end water; Lebrija's measured shrinkage curve reaches 1.95.
    The same cap converts the Laikipia vertic clays (prepare_willard.py);
    Vertisol_shrinkage/ has the memo and the cap test.
    The bulk-density column is the density at 33 kPa, as KSSL reports its
    vertisols' density at 1/3 bar.
  * Texture: treatment means for 0-0.20 m from G. Martinez's thesis (Table
    2.2, Textura tomejil.xlsx; LT = CT, SD = DD): sand 8.30 / 8.01 %, clay
    55.20 / 55.99 %, silt by difference -> clay (sand quartiles 6.2-9.2 %,
    clear of the 4.5 % that would make a 55 % clay soil silty clay). Not
    measured per ring, so every ring of a treatment carries the same
    fractions.
  * DD11's 10 and 31.6 cm readings repeated DD10's in the file; the
    correct values (0.43047, 0.4009) come from an earlier version of it
    (K. Vanderlinden). Dropped points: the 4,700 cm reading shared exactly
    by DD10 and CT3 (no explanation found), and dry-end readings wetter
    than the sample's own 500 cm reading (7 points in 5 curves; impossible
    on a drying curve, WP4 readings below ~0.5 MPa are imprecise).
  * SOC 1.1 (DD) and 0.9 % (CT), 0-0.1 m (Martinez et al. 2012).

Donana (Donana/RetencionAbalarioAcebuche+Ks Philip-Dunne.xlsx; site notes
in LOCALIZACION Y DESCRIPCION PUNTOS.docx):
  * Volumetric water content, the mean of 2-3 repeated rings per depth (47
    rings), on sand boxes at 1-500 mbar (0.1-50 kPa); El Acebuche also has
    WP4 readings beyond. El Abalario: 11 depths, 20-200 cm and the bottom
    (> 200 cm, depth left empty); El Acebuche: 20-100 cm and the surface
    (5-10 cm, depth 7.5). profile_id is the site.
  * Bulk density of the 47 rings averages 1.70 g/cm3 (1.59-1.80); the
    rings are not mapped to depths, so every layer carries 1.70.
  * Particle size (GRANULOMETRIA ABALARIO-ACEBUCHE.xls, IGME Sevilla 2010;
    K. Vanderlinden, email 2026-10-21): cumulative % passing, sieves to
    0.063 mm and sedimentation below, at the same depths as the curves
    (El Abalario 20-200 cm and > 200 cm, El Acebuche 20-100 cm). USDA
    fractions: sand = 100 - passing 0.05 mm, clay = passing 0.002 mm. The
    El Acebuche surface curve (5-10 cm) has no analysis; it keeps
    texture_class "sand" ("practically all sand", K. Vanderlinden) with
    the fractions empty.
  * Ks: Philip-Dunne infiltrometer in the field at 5-10 cm, repeated at
    different points, not on the retention samples (El Abalario 7 readings,
    El Acebuche 4). Kept as each site's median in ksat_pd_site_cmh for
    information; ksat_cmh stays empty.
  * Coordinates: UTM zone 29, taken as ED50 (X_29, Y_29 in the site notes).
  * No organic matter data.

Lebrija (Lebrija/CR_Lebrija 1054A_RAfa.xlsx; K. Vanderlinden, email
2026-10-09):
  * 8 sampling points in one field of Sector B-XII, lower Guadalquivir
    (Marisma recuperada, clay). One curve per point, so profile_id is the
    point; 12-14 readings each, -1 to -2.5e6 cm.
  * Gravimetric water content (theta_g) plus the volumetric water content
    (theta_v) the author derived with a shrinkage curve (Cornelis et al.
    2006) fitted to field void ratios from Andreu Caceres (1992, PhD thesis,
    Univ. Cordoba, Fig. 21; K. Vanderlinden, email 2026-10-21). theta_v is used as
    given; the bulk density it implies, theta_v/theta_g, runs 1.02 (wet) to
    1.95 g/cm3 and is capped at 1.951 at the dry end in every point, i.e.
    the fitted shrinkage limit. Unlike Tomejil, this is a measured
    shrinkage curve rather than the normal-shrinkage model: its specific
    volume sits 4-32 % above the normal-shrinkage line 1/2.65 + w, so
    structural and residual shrinkage are included. The bulk-density column
    is that implied density at 33 kPa, as for Tomejil.
  * Not a check: the wettest theta_v is 0.94 times the porosity of its own
    implied wet density in every point, but only because the fit, which
    stops at about 0.38 g/g, is extended along a normal-shrinkage line that
    holds saturation at 0.94 (data/andalusia_raw/Lebrija/NOTES.md).
  * Against the thesis (pressure plates, 15-60 cm): water near 33 kPa agrees
    (about 0.33 g/g in both); at 1500 kPa the thesis has about 0.32 g/g
    against 0.16-0.18 here (WP4C dew point). The gap is in gravimetric
    water, so not a matter of density. The thesis Ks, 0.46 cm/h on
    undisturbed columns at 5-40 cm, is not in the table.
  * Texture, organic matter, CEC, pH and saturation-extract EC per point
    from "Datos analiticos" (0-30 cm). Clay 56.7-68.8 %, sand 6.5-9.9 % ->
    clay in all of them. Point 5 was not analysed, so it is kept out of the
    table; its curve is still read by measured_points for prediction.
  * depth_cm 15, the midpoint of the 0-30 cm sampled layer.
  * No Ks. Coordinates: one pair for the field, UTM zone 29, taken as ED50
    as for Donana (36.98 N, 6.04 W).

Curves are fitted with the tool's own routine on the points up to 1500 kPa,
as for Arizona (prepare_babaeian_az.py); measured_points(site, None) keeps
the whole dry end. Tested both ways (new source, default arm): Tomejil with
its WP4 readings beyond 1500 kPa falls from 87 % to 0 % exact (every ring
read as silty clay loam); Setenil, where only 6 of 44 rings have a WP4
reading below 1500 kPa, scores the same exact class either way (31.8 %) and
a better group without the rest (54.5 vs 43.2 %).
"""

import os
import sys

import numpy as np
import pandas as pd

import free_m
import swcc_texture as st

RAW = os.path.join("data", "andalusia_raw")
SETENIL_XLSX = os.path.join(RAW, "Setenil", "SWRC_Setenil_Rafa.xlsx")
SETENIL_XY = os.path.join(RAW, "Setenil", "Puntos muestreo suelo_X_C.xls")
TOMEJIL_XLS = os.path.join(RAW, "Tomejil", "SWRCdata Tomejil.xls")
TOMEJIL_XY = os.path.join(RAW, "Tomejil", "coordenadas.xlsx")
DONANA_XLSX = os.path.join(RAW, "Do\u00f1ana", "RetencionAbalarioAcebuche+Ks Philip-Dunne.xlsx")
DONANA_PSD = os.path.join(RAW, "Do\u00f1ana", "GRANULOMETR\u00cdA ABALARIO-ACEBUCHE.xls")
LEBRIJA_XLSX = os.path.join(RAW, "Lebrija", "CR_Lebrija 1054A_RAfa.xlsx")
OUT = {"setenil": os.path.join("data", "setenil_reference.csv"),
       "tomejil": os.path.join("data", "tomejil_reference.csv"),
       "donana": os.path.join("data", "donana_reference.csv"),
       "lebrija": os.path.join("data", "lebrija_reference.csv")}
MBAR_TO_KPA = 0.1
CM_TO_KPA = 0.0980665
OM_PER_OC = 1.724
MIN_POINTS = 5
MAX_RMSE = 0.03
MAX_KPA = 1500.0
TOMEJIL_BD_MAX = 1.75          # g/cm3, median oven-dry clod density of KSSL Vertisols (see Tomejil notes)
TOMEJIL_SAND = {"DD": 8.01, "CT": 8.30}
TOMEJIL_CLAY = {"DD": 55.99, "CT": 55.20}
TOMEJIL_OC = {"DD": 1.1, "CT": 0.9}
# (sample, h cm) -> corrected reading, from an earlier version of the file.
TOMEJIL_FIXED = {("DD11", 10.0): 0.43047, ("DD11", 31.6): 0.4009}
# (sample, h cm) readings dropped: identical in two samples, unexplained.
TOMEJIL_COPIES = {("DD10", 4700.0), ("CT3", 4700.0)}
DONANA_BD = 1.70
DONANA_SITES = {"AB": ("El Abalario", 705321, 4109711),
                "AC": ("El Acebuche", 716291, 4102657)}
LEBRIJA_XY = (763950, 4097476)   # UTM zone 29, one pair for the field
LEBRIJA_DEPTH = 15.0             # midpoint of the sampled 0-30 cm layer


def _setenil_curves():
    d = pd.read_excel(SETENIL_XLSX, sheet_name="Final", header=None)
    out = {}
    for ring in range(1, 45):
        c = 2 * (ring - 1)
        assert str(d.iloc[0, c + 1]) == f"peso_{ring}"
        h = pd.to_numeric(d.iloc[1:, c], errors="coerce").to_numpy(float)
        w = pd.to_numeric(d.iloc[1:, c + 1], errors="coerce").to_numpy(float)
        ok = np.isfinite(h) & np.isfinite(w)
        out[ring] = (h[ok], w[ok])
    return out


def _setenil_rings():
    a = pd.read_excel(SETENIL_XLSX, sheet_name="Datos Agrario_anillos")
    b = pd.read_excel(SETENIL_XLSX, sheet_name="Datos Agrario bolsas")
    k = pd.read_excel(SETENIL_XLSX, sheet_name="CR_columna")
    r = pd.DataFrame({
        "ring": a.iloc[:, 0].astype(int), "point": a.iloc[:, 1].astype(int),
        "row": a.iloc[:, 2].astype(int), "depth": a.iloc[:, 3].astype(float),
        "clay": a["Arcilla"], "sand": a["Arena"], "silt": a["Limo"],
        "bd": a["densidad"], "stones": a["%piedras (w)"],
        "om": b.iloc[:, 15].astype(float), "ks": k["Ks(cm/h)"]})
    assert (b.iloc[:, 0].astype(int) == r.ring).all()
    assert (k.iloc[:, 0].astype(int) == r.ring).all()
    return r


def _tomejil_curves():
    x = pd.ExcelFile(TOMEJIL_XLS)
    out = {}
    for treat in ("DD", "CT"):
        d = x.parse(treat, header=None)
        for c in range(0, d.shape[1], 2):
            num = int(str(d.iloc[0, c + 1]).replace("humedad", ""))
            name = f"{treat}{num}"
            h = pd.to_numeric(d.iloc[1:, c], errors="coerce").to_numpy(float)
            w = pd.to_numeric(d.iloc[1:, c + 1], errors="coerce").to_numpy(float)
            ok = np.isfinite(h) & np.isfinite(w)
            ok &= np.array([(name, v) not in TOMEJIL_COPIES for v in h])
            h, w = h[ok], w[ok]
            w = np.array([TOMEJIL_FIXED.get((name, v), x) for v, x in zip(h, w)])
            w500 = w[np.argmin(np.abs(h - 500.0))]
            keep = (h <= 500.0) | (w <= w500)
            out[name] = (treat, num, h[keep], w[keep], int((~keep).sum()))
    return out


def shrinkage_bd(w, bd_max=TOMEJIL_BD_MAX, rho_s=2.65):
    """Bulk density at gravimetric water content w under normal shrinkage,
    capped at bd_max (the shrinkage limit)."""
    return np.minimum(bd_max, 1.0 / (1.0 / rho_s + np.asarray(w, float)))


def _donana_curves():
    """{(site, depth label): (depth_cm, h kPa, theta)}."""
    x = pd.ExcelFile(DONANA_XLSX)
    out = {}
    d = x.parse("ABALARIO", header=None)
    hdr = d.iloc[13].tolist()
    assert str(hdr[0]).startswith("Presi")
    h = pd.to_numeric(d.iloc[14:, 0], errors="coerce").to_numpy(float)
    for c in range(1, 12):
        lab = str(hdr[c]).strip()
        th = pd.to_numeric(d.iloc[14:, c], errors="coerce").to_numpy(float)
        ok = np.isfinite(h) & np.isfinite(th)
        depth = np.nan if lab.startswith("Fondo") else float(lab.split()[0])
        key = "FONDO" if lab.startswith("Fondo") else f"{int(depth):03d}"
        out[("AB", key)] = (depth, h[ok] * MBAR_TO_KPA, th[ok])
    d = x.parse("ACEBUCHE", header=None)
    hdr = d.iloc[13].tolist()
    for c in range(0, 12, 2):
        assert str(hdr[c]).startswith("PRESI"), hdr[c]
        lab = str(hdr[c + 1]).strip()
        h = pd.to_numeric(d.iloc[14:, c], errors="coerce").to_numpy(float)
        th = pd.to_numeric(d.iloc[14:, c + 1], errors="coerce").to_numpy(float)
        ok = np.isfinite(h) & np.isfinite(th)
        h, th = h[ok], th[ok]
        th500 = th[np.argmin(np.abs(h - 500.0))]
        keep = (h <= 500.0) | (th <= th500)
        if "superficie" in lab:
            depth, key = 7.5, "SUP"
        else:
            depth = float(lab.split("-")[1].split()[0]); key = f"{int(depth):03d}"
        out[("AC", key)] = (depth, h[keep] * MBAR_TO_KPA, th[keep])
    return out


def _donana_psd():
    """{(site, key): (sand, silt, clay)} in %, USDA limits."""
    sh = pd.read_excel(DONANA_PSD, header=None)
    sizes = pd.to_numeric(sh.iloc[1:18, 0], errors="coerce").to_numpy(float)
    out = {}
    for site, cols in (("AC", range(1, 6)), ("AB", range(8, 19))):
        for c in cols:
            lab = str(sh.iloc[0, c]).strip()
            key = "FONDO" if lab.startswith(">") else f"{int(lab.split()[0]):03d}"
            passing = dict(zip(sizes, pd.to_numeric(sh.iloc[1:18, c]).to_numpy(float)))
            clay, fine = passing[0.002], passing[0.05]
            out[(site, key)] = (round(100.0 - fine, 2), round(fine - clay, 2), clay)
    return out


def _donana_ks():
    """Each site's Philip-Dunne readings, cm/h."""
    x = pd.ExcelFile(DONANA_XLSX)
    ks = {}
    for site, sheet, rows in (("AB", "ABALARIO", range(3, 10)),
                              ("AC", "ACEBUCHE", range(4, 9))):
        d = x.parse(sheet, header=None)
        v = pd.to_numeric(d.iloc[list(rows), 1], errors="coerce").to_numpy(float)
        ks[site] = v[np.isfinite(v)] * 3.6e5          # m/s -> cm/h
    return ks


def _lebrija_curves():
    """{point: (h cm, theta_g, theta_v)}; three columns per point."""
    d = pd.read_excel(LEBRIJA_XLSX, sheet_name="CR", header=None)
    out = {}
    for p in range(1, 9):
        c = 3 * (p - 1)
        assert str(d.iloc[1, c]) == f"h (cm){p}"
        v = [pd.to_numeric(d.iloc[2:, c + i], errors="coerce").to_numpy(float)
             for i in range(3)]
        ok = np.isfinite(v[0]) & np.isfinite(v[1]) & np.isfinite(v[2])
        out[p] = (np.abs(v[0][ok]), v[1][ok], v[2][ok])
    return out


def _lebrija_analysis():
    """Per point, for the 7 points that were analysed."""
    d = pd.read_excel(LEBRIJA_XLSX, sheet_name="Datos analíticos")
    d = d.dropna(subset=["ARCILLA"])
    return pd.DataFrame({
        "point": d["Punto"].astype(int), "sand": d["ARENA"].astype(float),
        "silt": d["LIMO"].astype(float), "clay": d["ARCILLA"].astype(float),
        "om": d["MO (%)"].astype(float)}).set_index("point")


def _utm30_ed50(x, y, zone=30):
    import rasterio.warp
    lon, lat = rasterio.warp.transform(f"EPSG:230{zone}", "EPSG:4326", list(x), list(y))
    return np.array(lat), np.array(lon)


def measured_points(site, max_kpa=MAX_KPA):
    """{layer_id: (h kPa, theta)} for every sample of a site, for tests that
    refit; max_kpa=None keeps the whole dry end."""
    out = {}
    if site == "setenil":
        rings = _setenil_rings().set_index("ring")
        for ring, (h, w) in _setenil_curves().items():
            out[f"SET_R{ring:02d}"] = (h * CM_TO_KPA, w * rings.bd[ring])
    elif site == "tomejil":
        for name, (treat, _, h, w, _) in _tomejil_curves().items():
            out[f"TOM_{name}"] = (h * CM_TO_KPA, w * shrinkage_bd(w))
    elif site == "donana":
        for (s_, key), (_, h, th) in _donana_curves().items():
            out[f"DON_{s_}_{key}"] = (h, th)
    elif site == "lebrija":
        for p, (h, _, tv) in _lebrija_curves().items():
            out[f"LEB_P{p}"] = (h * CM_TO_KPA, tv)
    else:
        raise ValueError(site)
    if max_kpa is not None:
        out = {k: (h[h <= max_kpa], th[h <= max_kpa]) for k, (h, th) in out.items()}
    return out


def fit(h, theta):
    popt, _ = st.fit_vg(h, theta)
    tr, ts, la, ln1 = popt
    alpha, n = 10.0 ** la, 1.0 + 10.0 ** ln1
    rmse = float(np.sqrt(np.mean((st.vg_theta(h, tr, ts, alpha, n) - theta) ** 2)))
    fm = free_m.fit_free_m(h, theta,
                           fallback=free_m.mualem_fallback(tr, ts, alpha, n))
    return tr, ts, alpha, n, rmse, fm


def build_setenil():
    pts = measured_points("setenil")
    rings = _setenil_rings()
    xy = pd.read_excel(SETENIL_XY)
    lat, lon = _utm30_ed50(xy.iloc[:, 3], xy.iloc[:, 4])
    where = {(int(p), int(r)): (la, lo)
             for p, r, la, lo in zip(xy.iloc[:, 0], xy.iloc[:, 1], lat, lon)}
    rows, skip = [], dict(points=0, texture=0, fit=0, rmse=0)
    for r in rings.itertuples():
        lid = f"SET_R{r.ring:02d}"
        h, theta = pts[lid]
        if len(h) < MIN_POINTS:
            skip["points"] += 1
            continue
        tot = r.sand + r.silt + r.clay
        if not 95 <= tot <= 105:
            skip["texture"] += 1
            continue
        sand, silt, clay = (100 * v / tot for v in (r.sand, r.silt, r.clay))
        try:
            tr, ts, alpha, n, rmse, fm = fit(h, theta)
        except Exception:
            skip["fit"] += 1
            continue
        if not np.isfinite(rmse) or rmse > MAX_RMSE or n <= 1.0:
            skip["rmse"] += 1
            continue
        la, lo = where[(r.point, r.row)]
        rows.append(dict(
            layer_id=lid, profile_id=f"SET_P{r.point:02d}",
            texture_class=st.usda_class(sand, silt, clay),
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=sand, silt=silt, clay=clay, ksat_cmh=float(r.ks),
            depth_cm=float(r.depth), lat=la, lon=lo,
            source_db="Vanderlinden_Setenil", oc=r.om / OM_PER_OC, om=r.om,
            porosity=1 - r.bd / 2.65, bd=float(r.bd), rmse=rmse,
            n_points=len(h), **fm, sample_type="undisturbed",
            sample_type_source="intact rings: suction tables, WP4 and "
                               "laboratory permeameter (site notes)",
            position="inter-row" if r.row else "canopy",
            stones_pct=float(r.stones)))
    return pd.DataFrame(rows), skip


def build_tomejil():
    pts = measured_points("tomejil")
    curves = _tomejil_curves()
    xy = pd.read_excel(TOMEJIL_XY)
    lat, lon = _utm30_ed50(xy.iloc[:, 1], xy.iloc[:, 2])
    where = dict(zip(xy.iloc[:, 0].astype(int), zip(lat, lon)))
    rows, skip = [], dict(points=0, fit=0, rmse=0)
    for name, (treat, num, h_cm, w, dropped) in curves.items():
        lid = f"TOM_{name}"
        h, theta = pts[lid]
        w33 = np.interp(np.log10(330.0), np.log10(np.maximum(h_cm, 1e-3)), w)
        if len(h) < MIN_POINTS:
            skip["points"] += 1
            continue
        try:
            tr, ts, alpha, n, rmse, fm = fit(h, theta)
        except Exception:
            skip["fit"] += 1
            continue
        if not np.isfinite(rmse) or rmse > MAX_RMSE or n <= 1.0:
            skip["rmse"] += 1
            continue
        la, lo = where[num]
        bd, oc = float(shrinkage_bd(w33)), TOMEJIL_OC[treat]
        sand, clay = TOMEJIL_SAND[treat], TOMEJIL_CLAY[treat]
        rows.append(dict(
            layer_id=lid, profile_id=f"TOM_L{num:02d}",
            texture_class=st.usda_class(sand, 100 - sand - clay, clay),
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=sand, silt=100 - sand - clay, clay=clay,
            ksat_cmh=np.nan, depth_cm=2.5, lat=la, lon=lo,
            source_db="Vanderlinden_Tomejil", oc=oc, om=oc * OM_PER_OC,
            porosity=1 - bd / 2.65, bd=bd, rmse=rmse, n_points=len(h), **fm,
            sample_type="undisturbed",
            sample_type_source="intact 5 cm rings (Vanderlinden et al. 2021)",
            treatment=treat, dry_end_dropped=dropped))
    return pd.DataFrame(rows), skip


def build_donana():
    pts = measured_points("donana")
    curves = _donana_curves()
    ks = _donana_ks()
    psd = _donana_psd()
    rows, skip = [], dict(points=0, fit=0, rmse=0)
    for (site, key), (depth, _, _) in curves.items():
        lid = f"DON_{site}_{key}"
        h, theta = pts[lid]
        if len(h) < MIN_POINTS:
            skip["points"] += 1
            continue
        try:
            tr, ts, alpha, n, rmse, fm = fit(h, theta)
        except Exception:
            skip["fit"] += 1
            continue
        if not np.isfinite(rmse) or rmse > MAX_RMSE or n <= 1.0:
            skip["rmse"] += 1
            continue
        name, x, y = DONANA_SITES[site]
        (la,), (lo,) = _utm30_ed50([x], [y], zone=29)
        sand, silt, clay = psd.get((site, key), (np.nan,) * 3)
        rows.append(dict(
            layer_id=lid, profile_id=f"DON_{site}",
            texture_class=st.usda_class(sand, silt, clay) if np.isfinite(sand) else "sand",
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=sand, silt=silt, clay=clay, ksat_cmh=np.nan,
            depth_cm=depth, lat=la, lon=lo, source_db="Vanderlinden_Donana",
            oc=np.nan, om=np.nan, porosity=1 - DONANA_BD / 2.65, bd=DONANA_BD,
            rmse=rmse, n_points=len(h), **fm, sample_type="undisturbed",
            sample_type_source="rings on sand boxes (+ WP4 at El Acebuche); "
                               "mean of 2-3 rings per depth",
            site=name, ksat_pd_site_cmh=float(np.median(ks[site])),
            ksat_pd_n=len(ks[site])))
    return pd.DataFrame(rows), skip


def build_lebrija():
    pts = measured_points("lebrija")
    curves = _lebrija_curves()
    analysis = _lebrija_analysis()
    (la,), (lo,) = _utm30_ed50([LEBRIJA_XY[0]], [LEBRIJA_XY[1]], zone=29)
    rows, skip = [], dict(points=0, fit=0, rmse=0, texture=0)
    for p, (h_cm, w, tv) in curves.items():
        lid = f"LEB_P{p}"
        h, theta = pts[lid]
        if p not in analysis.index:
            skip["texture"] += 1
            continue
        if len(h) < MIN_POINTS:
            skip["points"] += 1
            continue
        try:
            tr, ts, alpha, n, rmse, fm = fit(h, theta)
        except Exception:
            skip["fit"] += 1
            continue
        if not np.isfinite(rmse) or rmse > MAX_RMSE or n <= 1.0:
            skip["rmse"] += 1
            continue
        # Implied density at 33 kPa, as for Tomejil.
        bd = float(np.interp(np.log10(330.0),
                             np.log10(np.maximum(h_cm, 1e-3)), tv / w))
        a = analysis.loc[p]
        rows.append(dict(
            layer_id=lid, profile_id=lid,
            texture_class=st.usda_class(a.sand, a.silt, a.clay),
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=a.sand, silt=a.silt, clay=a.clay,
            ksat_cmh=np.nan, depth_cm=LEBRIJA_DEPTH, lat=la, lon=lo,
            source_db="Vanderlinden_Lebrija", oc=a.om / OM_PER_OC, om=a.om,
            porosity=1 - bd / 2.65, bd=bd, rmse=rmse, n_points=len(h), **fm,
            sample_type="undisturbed",
            sample_type_source="theta_v from a shrinkage curve fitted to "
                               "void ratios (Andreu Caceres 1992)"))
    return pd.DataFrame(rows), skip


def main():
    if not os.path.isdir(RAW):
        raise SystemExit(f"{RAW} not found")
    sites = sys.argv[1:] or ["setenil", "tomejil", "donana", "lebrija"]
    for site in sites:
        out, skip = {"setenil": build_setenil, "tomejil": build_tomejil,
                     "donana": build_donana, "lebrija": build_lebrija}[site]()
        out.to_csv(OUT[site], index=False)
        print(f"{site}: {len(out)} layers -> {OUT[site]}   (skipped: {skip})")
        print(f"  median fit RMSE {out.rmse.median():.4f} (max {out.rmse.max():.4f}),"
              f" points per curve {out.n_points.min()}-{out.n_points.max()}")
        if out.ksat_cmh.notna().any():
            print(f"  Ks median {out.ksat_cmh.median():.2f} cm/h "
                  f"(range {out.ksat_cmh.min():.2f}-{out.ksat_cmh.max():.1f})")
        print("  " + out.texture_class.value_counts().to_string().replace("\n", "\n  "))


if __name__ == "__main__":
    main()
