"""One-time preprocessing of the Canary Islands andic soils (Armas Espinel).

Armas Espinel, S. (2013). Contribución al estudio de las propiedades físicas
de suelos ándicos de las Islas Canarias. PhD thesis, Universidad de La Laguna
(directors J. M. Hernández Moreno and C. M. Regalado), Serie Tesis Doctorales,
Ciencias y Tecnologías 42. Summarised in Armas-Espinel, S., Hernández-Moreno,
J. M., Muñoz-Carpena, R. and Regalado, C. M. (2003). Physical properties of
"sorriba"-cultivated volcanic soils from Tenerife in relation to andic
diagnostic parameters. Geoderma 117:297-311.

Reads data/armas_raw/cp471.pdf (the thesis) and writes
data/armas_reference.csv and data/armas_points.csv (the measured retention
points). Needs poppler (pdftotext, pdftocairo). All three are git-ignored.

  * Retention: undisturbed rings (96.6 cm3), saturated from below; Tempe
    cells from 0.2 to 90 kPa and Richards plates at 100, 500 and 1500 kPa
    (thesis IV.3.2). The measured points are not tabulated, but Anexo 3,
    Figura 1 (PDF pages 375-385) draws them as vector graphics, so they are
    read exactly: every grey marker is located and converted to pF and
    volumetric water content with each panel's own tick labels. pF below
    0.1 is taken as saturation (h = 0). Checked against the thesis: water
    content at 1500 kPa matches Anexo 1 Tablas 8-9 (median difference -0.3
    points), and refitting reproduces the thesis's own van Genuchten fits
    (Anexo 3 Tabla 5, m = 1-1/n: median alpha ratio 1.03, n-1 ratio 0.97).
  * Texture: Bouyoucos hydrometer after hexametaphosphate (HMP) dispersion,
    sand sieved at 0.2 and 0.05 mm, USDA limits (Anexo 1 Tablas 8-9); used as
    texture_class, taking another sampling date of the same sample where the
    curve's own date lacks it. Resin dispersion, which disperses andic
    material far more completely (Anexo 2 Tabla 1), is kept as
    sand_resin / silt_resin / clay_resin.
  * Ks: constant-head permeameter on the undisturbed rings (KsL, mm/h ->
    cm/h), same sampling date only, since it varied between campaigns.
  * Bulk density: cylinder method (Tablas 8-9), same date. Organic matter
    (Tablas 6-7, g/kg) -> om %, oc = om / 1.724.
  * Depth: mid-depth of the sampled layer (Tablas 4-5).
  * Coordinates: approximate, from the locality of each soil (Tabla 3); the
    thesis gives none.
"""

import html
import os
import re
import subprocess
import tempfile

import numpy as np
import pandas as pd

import free_m
import swcc_texture as st

PDF = os.path.join("data", "armas_raw", "cp471.pdf")
OUT = os.path.join("data", "armas_reference.csv")
POINTS = os.path.join("data", "armas_points.csv")
CM_TO_KPA = 0.0980665
FIG_PAGES = range(375, 386)          # Anexo 3, Figura 1
OM_PER_OC = 1.724
MIN_POINTS = 5
MAX_RMSE = 0.03

# Approximate locality coordinates (Tabla 3 gives the locality, not a point).
LOCALITY = {
    "N1": (28.83, -17.80, "La Palma", "Laguna de Barlovento"),
    "N2": (28.47, -16.40, "Tenerife", "Aguagarcía, Tacoronte"),
    "N3": (28.53, -16.28, "Tenerife", "Las Mercedes, La Laguna"),
    "N4": (28.43, -16.38, "Tenerife", "La Esperanza, El Rosario"),
    "C1": (28.52, -16.38, "Tenerife", "Valle Guerra, La Laguna"),
    "C2": (28.08, -16.65, "Tenerife", "Valle de San Lorenzo, Arona"),
    "C3": (28.08, -16.65, "Tenerife", "Valle de San Lorenzo, Arona"),
    "C4": (28.08, -16.65, "Tenerife", "Valle de San Lorenzo, Arona"),
    "C10": (28.52, -16.38, "Tenerife", "Valle Guerra, La Laguna"),
    "C11": (28.13, -15.63, "Gran Canaria", "Santa María de Guía"),
    "C12": (28.12, -15.52, "Gran Canaria", "Arucas"),
}

DATE = re.compile(r"^(Ene|Feb|Marz?|Abr|May|Jun|Jul|Ago|Sept?|Oct|Nov|Dic)\.?(\d{2})$", re.I)
SAMPLE = re.compile(r"^[CN]\d+[SP]\d?-\d+['’´]?$")


def norm_date(t):
    m = DATE.match(t)
    mon = m.group(1).lower().replace("marz", "mar").replace("sept", "sep")
    return f"{mon}{m.group(2)}"


def norm_sample(t):
    return re.sub(r"[’´]", "'", html.unescape(t))


def num(t):
    if t.strip() in ("n.d.", "nd", "n.d", "-", "–"):
        return np.nan
    if re.fullmatch(r"-?[\d.]*\d,\d+|-?\d+", t.strip()):
        return float(t.replace(".", "").replace(",", "."))
    return None


# --- measured points from the vector figures ---------------------------------

def _page_words(page, tmp):
    out = os.path.join(tmp, f"w{page}.html")
    subprocess.run(["pdftotext", "-f", str(page), "-l", str(page), "-bbox",
                    PDF, out], check=True)
    w = re.findall(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" '
                   r'yMax="([\d.]+)">([^<]*)</word>', open(out, encoding="utf-8").read())
    w = pd.DataFrame([(float(a), float(b), float(c), float(d), html.unescape(t))
                      for a, b, c, d, t in w], columns=["x0", "y0", "x1", "y1", "t"])
    w["v"] = w.t.map(lambda t: num(t) if re.fullmatch(r"\d+", t) else None)
    return w


def _page_markers(page, tmp):
    """Centres of the grey 'medido' markers. Each is drawn as a run of small
    triangles; a run ends when adding a triangle would exceed marker size."""
    out = os.path.join(tmp, f"p{page}.svg")
    subprocess.run(["pdftocairo", "-svg", "-f", str(page), "-l", str(page),
                    PDF, out], check=True)
    tris = re.findall(r'<path fill-rule="nonzero" fill="rgb\(50%, 50%, 50%\)"'
                      r'[^>]*d="([^"]*)"', open(out).read())
    groups, cur = [], []
    for d in tris:
        p = np.array(re.findall(r"([\d.]+) ([\d.]+)", d), float)
        if cur:
            a = np.vstack(cur + [p])
            if np.ptp(a[:, 0]) > 5.0 or np.ptp(a[:, 1]) > 5.0:
                groups.append(np.vstack(cur))
                cur = []
        cur.append(p)
    if cur:
        groups.append(np.vstack(cur))
    return np.array([[(g[:, 0].min() + g[:, 0].max()) / 2,
                      (g[:, 1].min() + g[:, 1].max()) / 2] for g in groups])


def measured_points():
    """{layer_id: (h kPa, theta)} for every curve, from data/armas_points.csv."""
    pts = pd.read_csv(POINTS)
    return {f"AR_{s}_{d}": (np.where(g.pF < 0.1, 0.0, 10.0 ** g.pF * CM_TO_KPA),
                            g.theta_pct.to_numpy() / 100.0)
            for (s, d), g in pts.groupby(["sample", "date"])}


def extract_points():
    """DataFrame of (sample, date, pF, theta_pct) read from Figura 1."""
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for page in FIG_PAGES:
            w, mk = _page_words(page, tmp), _page_markers(page, tmp)
            for _, r in w[w.t == "MUESTRA:"].iterrows():
                line = w[(abs(w.y0 - r.y0) < 2.5) & (w.x0 > r.x1)].sort_values("x0")
                sample = norm_sample(line.iloc[0].t)
                fecha = line[line.t == "FECHA:"]
                date = (line[line.x0 > fecha.x1.iloc[0]].iloc[0].t.lower()
                        .replace("sept", "sep").replace("marz", "mar")
                        if len(fecha) else None)
                yl = w[w.v.notna() & (w.x1 < r.x0 + 1) & (w.x1 > r.x0 - 25)
                       & (w.y0 > r.y1) & (w.y0 < r.y1 + 150)].drop_duplicates("v")
                xl = w[w.v.isin([0, 1, 2, 3, 4, 5]) & (w.x0 > r.x0 - 8)
                       & (w.x0 < r.x0 + 150) & (w.y0 > r.y1 + 100)
                       & (w.y0 < r.y1 + 160)]
                if len(xl):
                    xl = xl[abs(xl.y0 - xl.y0.max()) < 2].drop_duplicates("v")
                if len(yl) < 3 or len(xl) < 4:
                    raise RuntimeError(f"axis calibration failed: page {page}, {sample}")
                xc, yc = (xl.x0 + xl.x1) / 2, (yl.y0 + yl.y1) / 2
                ax, bx = np.polyfit(xc, xl.v.astype(float), 1)
                ay, by = np.polyfit(yc, yl.v.astype(float), 1)
                inside = mk[(mk[:, 0] > xc.min() - 3) & (mk[:, 0] < xc.max() + 1)
                            & (mk[:, 1] > yc.min() - 3) & (mk[:, 1] < yc.max() + 3)]
                for x, y in inside:
                    rows.append(dict(page=page, sample=sample, date=date,
                                     pF=ax * x + bx, theta_pct=ay * y + by))
    pts = pd.DataFrame(rows)
    # One curve (N1P2-1, Sept. 2004) is drawn twice; keep its first figure.
    first = pts.groupby(["sample", "date"]).page.transform("min")
    return pts[pts.page == first].drop(columns="page").reset_index(drop=True)


# --- tables ------------------------------------------------------------------

def _table_rows(pages, page_range):
    for i in page_range:
        for line in pages[i - 1].splitlines():
            tok = line.split()
            si = next((k for k, x in enumerate(tok) if SAMPLE.match(x)), None)
            if si is None:
                continue
            di = next((k for k, x in enumerate(tok[:si]) if DATE.match(x)), None)
            yield (norm_date(tok[di]) if di is not None else None,
                   norm_sample(tok[si]), tok[si + 1:])


def tables():
    txt = subprocess.run(["pdftotext", "-layout", PDF, "-"], check=True,
                         capture_output=True).stdout.decode("utf-8", "ignore")
    pages = txt.split("\f")
    depth = {}
    for i in range(31, 36):                                # Tablas 4-5
        for line in pages[i - 1].splitlines():
            tok = line.split()
            if len(tok) >= 3 and SAMPLE.match(tok[-1]):
                m = re.fullmatch(r"(\d+)-(\d+)", tok[-3]) or re.fullmatch(r"(\d+)", tok[-3])
                if m:
                    lo = float(m.group(1))
                    hi = float(m.group(2)) if m.lastindex == 2 else lo
                    depth[norm_sample(tok[-1])] = (lo + hi) / 2
    phys = []                                              # Anexo 1 Tablas 8-9
    for date, s, a in _table_rows(pages, range(334, 349)):
        v = [num(x) for x in a]
        if len(a) >= 10 and all(x is not None for x in v[:7]):
            phys.append(dict(date=date, sample=s, bd=v[0], th1500=v[3], clay=v[4],
                             silt=v[5], sand=v[6],
                             ksl=v[8] if v[8] is not None else np.nan))
    chem = []                                              # Anexo 1 Tablas 6-7
    for date, s, a in _table_rows(pages, range(328, 334)):
        v = [num(x) for x in a]
        if len(v) >= 12 and all(x is not None for x in v[:12]):
            chem.append(dict(date=date, sample=s, om_gkg=v[9]))
    disp = []                                              # Anexo 2 Tabla 1
    for date, s, a in _table_rows(pages, range(349, 350)):
        v = [num(x) for x in a]
        if len(a) >= 8 and all(x is not None for x in v[:3] + v[4:7]):
            disp.append(dict(sample=s, clay_resin=v[4], silt_resin=v[5], sand_resin=v[6]))
    return depth, pd.DataFrame(phys), pd.DataFrame(chem), pd.DataFrame(disp)


def main():
    if not os.path.exists(PDF):
        raise SystemExit(f"{PDF} not found")
    pts = extract_points()
    pts.to_csv(POINTS, index=False)
    depth, phys, chem, disp = tables()
    base = lambda s: s.replace("'", "")

    def pick(df, s, date, col, any_date=True):
        d = df[df["sample"].map(base) == base(s)]
        same = d[(d.date == date) & d[col].notna()]
        if len(same):
            return same.iloc[0]
        d = d[d[col].notna()]
        return d.iloc[0] if (any_date and len(d)) else None

    rows, skip = [], dict(points=0, texture=0, fit=0, rmse=0)
    for (s, date), g in pts.groupby(["sample", "date"]):
        h = np.where(g.pF < 0.1, 0.0, 10.0 ** g.pF * CM_TO_KPA)
        theta = g.theta_pct.to_numpy() / 100.0
        if len(h) < MIN_POINTS:
            skip["points"] += 1
            continue
        tex = pick(phys, s, date, "clay")
        if tex is None:
            skip["texture"] += 1
            continue
        tot = tex.sand + tex.silt + tex.clay
        sand, silt, clay = (100 * tex.sand / tot, 100 * tex.silt / tot,
                            100 * tex.clay / tot)
        try:
            popt, _ = st.fit_vg(h, theta)
        except Exception:
            skip["fit"] += 1
            continue
        tr, ts, la, ln1 = popt
        alpha, n = 10.0 ** la, 1.0 + 10.0 ** ln1
        rmse = float(np.sqrt(np.mean((st.vg_theta(h, tr, ts, alpha, n) - theta) ** 2)))
        if not np.isfinite(rmse) or rmse > MAX_RMSE or n <= 1.0:
            skip["rmse"] += 1
            continue
        fm = free_m.fit_free_m(h, theta,
                               fallback=free_m.mualem_fallback(tr, ts, alpha, n))
        soil = re.match(r"([CN]\d+)", s).group(1)
        lat, lon, island, locality = LOCALITY.get(soil, (np.nan, np.nan, "", ""))
        phy = pick(phys, s, date, "bd", any_date=False)
        ks = pick(phys, s, date, "ksl", any_date=False)
        om = pick(chem, s, date, "om_gkg")
        rz = disp[disp["sample"].map(base) == base(s)]
        resin = {}
        if len(rz):
            r = rz.iloc[0]
            rt = r.sand_resin + r.silt_resin + r.clay_resin
            resin = {c: 100 * r[c] / rt for c in ("sand_resin", "silt_resin", "clay_resin")}
        rows.append(dict(
            layer_id=f"AR_{s}_{date}",
            # Samples of one soil come from one field (or profile), so the
            # field is the site for grouped hold-outs.
            profile_id=f"AR_{soil}",
            texture_class=st.usda_class(sand, silt, clay),
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=sand, silt=silt, clay=clay,
            ksat_cmh=ks.ksl / 10.0 if ks is not None else np.nan,
            depth_cm=depth.get(base(s), np.nan),
            lat=lat, lon=lon, source_db="Armas_Canarias",
            oc=om.om_gkg / 10.0 / OM_PER_OC if om is not None else np.nan,
            om=om.om_gkg / 10.0 if om is not None else np.nan,
            porosity=np.nan, bd=phy.bd if phy is not None else np.nan,
            rmse=rmse, n_points=len(h), **fm,
            sample_type="undisturbed",
            sample_type_source="Armas Espinel (2013): undisturbed rings, Tempe "
                               "cells and pressure plates; Ks on the rings",
            island=island, locality=locality,
            sand_resin=resin.get("sand_resin", np.nan),
            silt_resin=resin.get("silt_resin", np.nan),
            clay_resin=resin.get("clay_resin", np.nan)))

    out = pd.DataFrame(rows)
    out["texture_class_resin"] = [
        st.usda_class(a, b, c) if np.isfinite(c) else None
        for a, b, c in out[["sand_resin", "silt_resin", "clay_resin"]].to_numpy()]
    out.to_csv(OUT, index=False)
    print(f"curves read:      {pts.groupby(['sample', 'date']).ngroups} "
          f"({len(pts)} points)")
    print(f"reference layers: {len(out)}   (skipped: {skip})")
    print(f"  soils/fields:   {out.profile_id.nunique()}  "
          f"{out.profile_id.str[3:].value_counts().to_dict()}")
    print(f"  with ksat:      {out.ksat_cmh.notna().sum()} (median "
          f"{out.ksat_cmh.median():.1f} cm/h)")
    print(f"  with bd:        {out.bd.notna().sum()} (median {out.bd.median():.2f})")
    print(f"  with om:        {out.om.notna().sum()} (median {out.om.median():.1f} %)")
    print(f"  with resin texture: {out.clay_resin.notna().sum()}")
    print(f"  median fit RMSE: {out.rmse.median():.4f}")
    print("  class (HMP):  ", out.texture_class.value_counts().to_dict())
    print("  class (resin):", out.texture_class_resin.value_counts().to_dict())


if __name__ == "__main__":
    main()
