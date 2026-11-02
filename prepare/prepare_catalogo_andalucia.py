"""Collate the Catalogo de suelos de Andalucia into one horizon table.

Junta de Andalucia, Consejeria de Medio Ambiente. Catalogo de suelos de
Andalucia (documentos tecnicos), one PDF per profile, public domain:
https://www.juntadeandalucia.es/medioambiente/documentos_tecnicos/catalogo2/catalogo_ind.html

Reads data/andalucia_catalogo_raw/*.pdf (62 profiles, almeria1 ... sevilla9,
plus metodologia.pdf) through `pdftotext -layout` and writes
data/catalogo_andalucia_horizons.csv: one row per horizon, the profile's
site description repeated on each row.

From metodologia.pdf:
  * Bulk density (DENS. APAR.) on undisturbed samples; total porosity (POR.)
    computed from bulk and particle density; Ks (COND. HIDR., cm/h) by water
    passed through saturated undisturbed cylinders.
  * Retention at three suctions only, weight %: 1/10 and 1/3 bar by
    centrifuge on sieved soil (the "moisture equivalent"), 15 bar on a
    Richards membrane, sieved soil.
  * Texture in four fractions, coarse sand 2-0.2, fine sand 0.2-0.05, silt
    0.05-0.002 and clay < 0.002 mm, mostly by hydrometer; "in some cases" the
    sand/silt limit was 0.02 mm instead of 0.05 mm, not said which.
  * Blank cells are determinations not made. Values are kept as printed:
    a blank silt beside sand + clay near 100 (H-04) is not filled in, and
    w10kpa = w33kpa in CA-05 AC is the sheet's.

Values are placed in their columns by the header's character positions in
the layout text (blank cells leave no token, so counting tokens fails).
Coordinates are UTM in km as printed; the zone is not stated. Western
profiles (X > 620 km, Huelva and western Seville/Cadiz) are taken as zone 29,
the rest as zone 30, both ED50.
"""

import glob
import os
import re
import subprocess

import numpy as np
import pandas as pd

RAW = os.path.join("data", "andalucia_catalogo_raw")
OUT = os.path.join("data", "catalogo_andalucia_horizons.csv")
PROVINCES = ["almeria", "cadiz", "cordoba", "granada", "huelva", "jaen",
             "malaga", "sevilla"]

# (header token, column name) per analytical table; the token's centre in the
# header line is the column's position.
T_CHEM = [("H20", "ph_h2o"), ("CLK", "ph_kcl"), ("C ORG.", "oc"),
          ("N TO.", "n_tot"), ("C/N", "cn"), ("P TOT.", "p_tot_ppm"),
          ("TOTAL", "fe_tot"), ("LIBRE", "fe_free"), ("AMORFO", "fe_amorph")]
T_CEC = [("(%)", "co3"), ("CA ", "ca"), ("MG ", "mg"), ("NA ", "na"),
         ("K ", "k"), ("H ", "h"), ("(SUMA)", "cec"), ("BAS.", "base_sat")]
T_PHYS = [("DENS.", "bd"), ("PORO", "porosity_pct"), ("COND.", "ksat_cmh"),
          ("1/10B", "w10kpa"), ("1/3B", "w33kpa"), ("15B", "w1500kpa"),
          ("2-0,2", "coarse_sand"), ("0,2-0,05", "fine_sand"),
          ("0,05-0,002", "silt"), ("< 0,002", "clay")]

NUM = re.compile(r"(?<![\w,])(\d+,\d+|\d+)(?![\w,])")
SITE_FIELDS = {
    "profile_id": r"PERFIL N\S*MERO:\s*(.+)",
    "land_use": r"USO ACTUAL:\s*(.+)",
    "elevation_m": r"ELEVACI\S*N:\s*(\d+)",
    "slope_pct": r"PENDIENTE:\s*([\d,.]+)",
    "erosion": r"EROSI\S*N:\s*(.+)",
    "drainage": r"DRENAJE:\s*(.+)",
    "parent_material": r"MATERIAL ORIGINAL:\s*(.+)",
    "position": r"POSICI\S*N FISIOGR\S*FICA:\s*(.+)",
    "soil_taxonomy": r"CLASIFICACI\S*N:\s*(.+)",
    "observer_date": r"OBSERVADOR Y FECHA:\s*(.+)",
}


def _num(s):
    return float(s.replace(",", "."))


def _text(pdf):
    return subprocess.run(["pdftotext", "-layout", pdf, "-"], check=True,
                          capture_output=True, text=True).stdout


def _centres(header_lines, spec):
    """Column centre per name, from the first header line holding its token."""
    out = {}
    for token, name in spec:
        for line in header_lines:
            i = line.find(token)
            if i >= 0:
                out[name] = i + len(token.strip()) / 2
                break
    return out


def _table(lines, start, n_header, spec, depth=False):
    """Rows of the table whose first header line is lines[start]:
    {horizon: {column: value}} in order, values placed by position."""
    cols = _centres(lines[start:start + n_header], spec)
    names, centres = list(cols), np.array(list(cols.values()))
    rows = []
    for line in lines[start + n_header:]:
        if not line.strip():
            if rows:
                break
            continue
        if not rows and ("(" in line or line.lstrip().startswith("-")):
            continue                      # units lines of the header
        if line.lstrip().startswith(("HORIZ", "NOTA")):
            break
        hz = line.split()[0] if line[0] != " " else ""
        rec = {"horizon": hz}
        body_start = len(hz)
        if depth:
            dm = re.match(r"\s*(\d+)(?![,\d])\s*-?\s*(\d+(?![,\d]))?\s*-?", line[len(hz):])
            if dm:
                rec["top_cm"] = float(dm.group(1))
                rec["bottom_cm"] = float(dm.group(2)) if dm.group(2) else np.nan
                body_start = len(hz) + dm.end()
        for t in NUM.finditer(line, body_start):
            if "," not in t.group(1) and depth:
                continue
            c = (t.start() + t.end()) / 2
            rec[names[int(np.argmin(abs(centres - c)))]] = _num(t.group(1))
        rows.append(rec)
    return rows


def parse_profile(pdf):
    text = _text(pdf)
    lines = text.splitlines()
    base = os.path.splitext(os.path.basename(pdf))[0]
    site = {"file": base, "province": re.sub(r"\d+", "", base),
            "soil_name": next(l.strip() for l in lines if l.strip())}
    for key, pat in SITE_FIELDS.items():
        m = re.search(pat, text)
        site[key] = m.group(1).strip() if m else None
    m = re.search(r"X\s*=\s*([\d.]+)\s*;\s*Y\s*=\s*([\d.]+)", text)
    site["utm_x_km"], site["utm_y_km"] = (float(m.group(1)), float(m.group(2))) if m else (np.nan, np.nan)
    for key in ("elevation_m", "slope_pct"):
        site[key] = _num(site[key].replace(".", ",")) if site[key] else np.nan
    m = re.search(r"(1[89]\d\d)", site["observer_date"] or "")
    site["year"] = int(m.group(1)) if m else np.nan

    i_chem = next(i for i, l in enumerate(lines) if l.lstrip().startswith("HORIZ.") and "PROF." in l)
    i_cec = next(i for i, l in enumerate(lines) if l.lstrip().startswith("HORIZ.") and "CO3" in l and "C.C.C." in l)
    i_phys = next(i for i, l in enumerate(lines) if l.lstrip().startswith("HORIZ.") and "DENS." in l)
    chem = _table(lines, i_chem, 2, T_CHEM, depth=True)
    cec = _table(lines, i_cec, 3, T_CEC)
    phys = _table(lines, i_phys, 3, T_PHYS)
    # Tables list the same horizons in the same order (one sheet prints no
    # horizon names in its physical table, so rows are matched by order).
    if not (len(chem) == len(cec) == len(phys)):
        raise ValueError(f"{base}: {len(chem)}/{len(cec)}/{len(phys)} horizon rows")
    rows = []
    for k, (a, b, c) in enumerate(zip(chem, cec, phys)):
        for other in (b, c):
            if other["horizon"] and other["horizon"] != a["horizon"]:
                raise ValueError(f"{base}: horizon {a['horizon']} vs {other['horizon']}")
        rec = {**site, "horizon_no": k + 1}
        rec.update(a)
        rec.update({key: v for key, v in b.items() if key != "horizon"})
        rec.update({key: v for key, v in c.items() if key != "horizon"})
        rows.append(rec)
    return rows


def main():
    pdfs = sorted(glob.glob(os.path.join(RAW, "*.pdf")))
    pdfs = [p for p in pdfs if re.fullmatch(r"[a-z]+\d+", os.path.splitext(os.path.basename(p))[0])]
    df = pd.DataFrame([r for p in pdfs for r in parse_profile(p)])

    for c in ("ph_h2o", "ph_kcl", "oc", "n_tot", "cn", "p_tot_ppm", "fe_tot", "fe_free",
              "fe_amorph", "co3", "ca", "mg", "na", "k", "h", "cec", "base_sat",
              "bd", "porosity_pct", "ksat_cmh", "w10kpa", "w33kpa", "w1500kpa",
              "coarse_sand", "fine_sand", "silt", "clay", "top_cm", "bottom_cm"):
        if c not in df:
            df[c] = np.nan
    df["sand"] = df.coarse_sand + df.fine_sand
    df["fraction_sum"] = df.sand + df.silt + df.clay
    df["depth_cm"] = np.where(df.bottom_cm.notna(), (df.top_cm + df.bottom_cm) / 2, df.top_cm)
    df["layer_id"] = ("CSA_" + df.profile_id.str.extract(r"^([A-Z]+-\d+)")[0]
                      + "_" + df.horizon_no.astype(str))
    df["n_retention"] = df[["w10kpa", "w33kpa", "w1500kpa"]].notna().sum(axis=1)

    import rasterio.warp
    df["utm_zone"] = np.where(df.utm_x_km > 620, 29, 30)
    df["lat"] = df["lon"] = np.nan
    for z, g in df.groupby("utm_zone"):
        lon, lat = rasterio.warp.transform(f"EPSG:230{z}", "EPSG:4326",
                                           list(g.utm_x_km * 1000), list(g.utm_y_km * 1000))
        df.loc[g.index, "lon"], df.loc[g.index, "lat"] = lon, lat
    # AL-03 prints Y = 4928.000, which is in France: a typo, left blank.
    off = ~df.lat.between(35.9, 38.8)
    df.loc[off, ["lat", "lon"]] = np.nan

    first = ["layer_id", "profile_id", "file", "province", "soil_name", "horizon_no",
             "horizon", "top_cm", "bottom_cm", "depth_cm"]
    df = df[first + [c for c in df.columns if c not in first]]
    df.to_csv(OUT, index=False)
    print(f"{len(pdfs)} profiles, {len(df)} horizons -> {OUT}")
    print(f"  with texture (4 fractions): {df.fraction_sum.notna().sum()}, "
          f"bulk density: {df.bd.notna().sum()}, Ks: {df.ksat_cmh.notna().sum()}, "
          f"porosity: {df.porosity_pct.notna().sum()}")
    print("  retention points per horizon: "
          + df.n_retention.value_counts().sort_index().to_string().replace("\n", ", "))


if __name__ == "__main__":
    main()
