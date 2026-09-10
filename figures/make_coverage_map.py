"""Coverage map of the reference tables, Equal Earth projection.

    .venv/bin/python figures/make_coverage_map.py                 # default reference
    .venv/bin/python figures/make_coverage_map.py public euhydi   # base + addition

With two arguments the base set is drawn in the accent colour and the added
set underneath it in a second colour, so what the addition contributes is
visible. Writes figures/fig2_coverage_map.svg (standalone, colours baked in)
and figures/fig2_coverage_map.fragment.svg (CSS variables, for workflow.html).

Only 2-degree cell counts are drawn, never individual sample locations, so the
figure can be shared even when a table is restricted.

Land outline: Natural Earth 110 m (public domain), downloaded once to
figures/ne_110m_land.geojson. Equal Earth: Savric, Patterson & Jenny (2018).
"""

import json
import math
import os
import sys
import urllib.request

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import swcc_texture as st                     # noqa: E402
from verify_region import EU_BBOX             # noqa: E402

LAND_URL = ("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/"
            "master/geojson/ne_110m_land.geojson")
LAND = os.path.join(HERE, "ne_110m_land.geojson")
BIN = 2.0
LAT_N, LAT_S = 84.0, -58.0
W, PAD_X, PAD_TOP = 1040.0, 14.0, 44.0
LABEL = {"euhydi": "EU-HYDI", "hohenbrink": "Hohenbrink", "kssl": "KSSL",
         "gshp": "GSHP", "merged": "Default reference",
         "public": "Public reference"}
BAKE = {"var(--accent)": "#0E6E73", "var(--add)": "#B8471E",
        "var(--muted)": "#5C6A70", "var(--sunken)": "#EAEEEC",
        "currentColor": "#141A1D"}
REGIONS = {"North America": (-170, -52, 15, 72),
           "Europe": (*EU_BBOX["lon"], *EU_BBOX["lat"]),
           "Asia": (32, 150, 5, 78), "South America": (-82, -34, -56, 13),
           "Africa": (-18, 52, -35, 37), "Australasia": (110, 180, -48, -8)}

A1, A2, A3, A4 = 1.340264, -0.081106, 0.000893, 0.003796
M = math.sqrt(3) / 2.0


def eqearth(lon, lat):
    lam = np.radians(np.asarray(lon, float))
    th = np.arcsin(np.clip(M * np.sin(np.radians(np.asarray(lat, float))), -1, 1))
    t2 = th * th
    t6 = t2 ** 3
    x = lam * np.cos(th) / (M * (A1 + 3 * A2 * t2 + t6 * (7 * A3 + 9 * A4 * t2)))
    return x, th * (A1 + A2 * t2 + t6 * (A3 + A4 * t2))


XW = eqearth([180.0], [0.0])[0][0]
SCALE = (W - 2 * PAD_X) / (2 * XW)
CY = PAD_TOP + eqearth([0.0], [90.0])[1][0] * SCALE


def px(lon, lat):
    x, y = eqearth(lon, lat)
    return W / 2 + x * SCALE, CY - y * SCALE


def simplify(pts, tol):
    """Douglas-Peucker on a list of (x, y)."""
    if len(pts) < 3:
        return pts
    a, b = np.array(pts[0]), np.array(pts[-1])
    ab = b - a
    q = np.array(pts) - a
    n = float(np.hypot(*ab))
    d = (np.hypot(q[:, 0], q[:, 1]) if n == 0
         else np.abs(ab[0] * q[:, 1] - ab[1] * q[:, 0]) / n)
    i = int(np.argmax(d))
    if d[i] > tol:
        return simplify(pts[:i + 1], tol)[:-1] + simplify(pts[i:], tol)
    return [pts[0], pts[-1]]


def path(xs, ys, close=False):
    return ("M" + "L".join(f"{x:.1f} {y:.1f}" for x, y in zip(xs, ys))
            + ("Z" if close else ""))


def land_paths():
    if not os.path.exists(LAND):
        urllib.request.urlretrieve(LAND_URL, LAND)
    out = []
    for f in json.load(open(LAND))["features"]:
        g = f["geometry"]
        rings = ([g["coordinates"][0]] if g["type"] == "Polygon"
                 else [p[0] for p in g["coordinates"]])
        for ring in rings:
            lat = [c[1] for c in ring]
            if max(lat) < LAT_S + 2:          # Antarctica: no samples
                continue
            xs, ys = px([c[0] for c in ring], lat)
            pts = simplify(list(zip(xs.tolist(), ys.tolist())), 0.9)
            if len(pts) >= 4:
                out.append(path(*zip(*pts), close=True))
    return out


def georef(df):
    d = df[df.lat.notna() & df.lon.notna()].copy()
    d["source_db"] = d.source_db.fillna("KSSL")
    return d


def cells(d):
    # Only cells inside the mapped latitudes: anything south of LAT_S would
    # otherwise plot below the map, in the legend.
    d = d[d.lat.between(LAT_S, LAT_N)]
    b = pd.DataFrame({"lon": (np.floor(d.lon / BIN) + .5) * BIN,
                      "lat": (np.floor(d.lat / BIN) + .5) * BIN})
    return b.groupby(["lon", "lat"]).size().reset_index(name="n")


def in_region(d, box):
    w, e, s, n = box
    return int((d.lon.between(w, e) & d.lat.between(s, n)).sum())


def main():
    base_name = sys.argv[1] if len(sys.argv) > 1 else "merged"
    add_name = sys.argv[2] if len(sys.argv) > 2 else None
    base = georef(st.load_reference_df(base_name))
    add = georef(st.load_reference_df(add_name)) if add_name else base.iloc[:0]
    allr = pd.concat([base, add], ignore_index=True)

    cb, ca = cells(base), cells(add)
    nmax = max(cb.n.max(), ca.n.max() if len(ca) else 0)
    rmax = 25.0

    def radius(n):
        return np.clip(rmax * np.sqrt(np.asarray(n, float) / nmax), 1.6, rmax)

    def circles(c):
        x, y = px(c.lon.to_numpy(), c.lat.to_numpy())
        return "\n".join(f'        <circle cx="{a:.1f}" cy="{b:.1f}" r="{r:.2f}"/>'
                         for a, b, r in zip(x, y, radius(c.n)))

    top = CY - eqearth([0.0], [LAT_N])[1][0] * SCALE
    bot = CY - eqearth([0.0], [LAT_S])[1][0] * SCALE
    vb_y, vb_h = top - 34, (bot - top) + 34 + 70

    grat = []
    for la in range(-40, 81, 40):
        lo = np.linspace(-180, 180, 90)
        grat.append(path(*px(lo, np.full_like(lo, la))))
    for lo in range(-120, 121, 60):
        la = np.linspace(LAT_S, LAT_N, 50)
        grat.append(path(*px(np.full_like(la, lo), la)))

    ew, ee = EU_BBOX["lon"]
    es, en = EU_BBOX["lat"]
    elon = np.concatenate([np.linspace(ew, ee, 40), np.full(24, ee),
                           np.linspace(ee, ew, 40), np.full(24, ew)])
    elat = np.concatenate([np.full(40, en), np.linspace(en, es, 24),
                           np.full(40, es), np.linspace(es, en, 24)])
    ex, ey = px(elon, elat)
    eu_b = in_region(base, REGIONS["Europe"])
    eu_t = in_region(allr, REGIONS["Europe"])
    eu_label = (f"Europe &#183; {eu_b:,} &#8594; {eu_t:,} layers" if add_name
                else f"Europe &#183; {eu_b:,} layers")
    eu_sub = f"{eu_t / len(allr) * 100:.1f} % of the reference"

    src = allr.source_db.value_counts()
    fl = int(src.get("Florida_database", 0))
    callouts = f'''
        <line x1="302" y1="192" x2="302" y2="228" stroke="currentColor" opacity=".55"/>
        <text x="296" y="243" text-anchor="end">Florida &#8212; {fl:,} layers,</text>
        <text x="296" y="257" text-anchor="end">{fl / len(allr) * 100:.0f} % of everything</text>
        <line x1="389" y1="336" x2="404" y2="356" stroke="currentColor" opacity=".55"/>
        <text x="408" y="360">Brazil (HYBRAS) {int(src.get("HYBRAS", 0)):,}</text>
        <line x1="619" y1="97" x2="660" y2="76" stroke="currentColor" opacity=".55"/>
        <text x="664" y="72">W. Russia {int(src.get("Russia_EGRPR", 0)):,}</text>
        <line x1="576" y1="298" x2="600" y2="322" stroke="currentColor" opacity=".55"/>
        <text x="604" y="326">Sub-Saharan Africa {int(src.get("AfSPDB", 0)):,}</text>
        <line x1="279" y1="149" x2="238" y2="120" stroke="currentColor" opacity=".55"/>
        <text x="234" y="116" text-anchor="end">KSSL, across the US {int(src.get("KSSL", 0)):,}</text>'''

    ly = bot + 34
    size_leg, lx = [], 20.0
    for cnt in (10, 100, 1000):
        r = float(radius([cnt])[0])
        size_leg.append((lx + r, ly, r, cnt))
        lx += 2 * r + 40
    size_svg = "\n".join(
        f'        <circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.2f}" fill="currentColor" '
        f'fill-opacity=".12" stroke="currentColor" stroke-opacity=".5"/>\n'
        f'        <text x="{x:.1f}" y="{y + 28:.0f}" text-anchor="middle">{c:,}</text>'
        for x, y, r, c in size_leg)
    sx = size_leg[-1][0] + size_leg[-1][2] + 14
    key_x = 640.0
    key = (f'        <circle cx="{key_x:.0f}" cy="{ly - 8:.0f}" r="6" fill="var(--accent)" '
           f'fill-opacity=".45" stroke="var(--accent)"/>\n'
           f'        <text x="{key_x + 14:.0f}" y="{ly - 4:.0f}">'
           f'{LABEL.get(base_name, base_name)} &#183; {len(base):,} georeferenced layers</text>')
    if add_name:
        key += (f'\n        <circle cx="{key_x:.0f}" cy="{ly + 14:.0f}" r="6" fill="var(--add)" '
                f'fill-opacity=".5" stroke="var(--add)"/>\n'
                f'        <text x="{key_x + 14:.0f}" y="{ly + 18:.0f}">'
                f'{LABEL.get(add_name, add_name)} &#183; +{len(add):,}</text>')

    aria = (f"Equal Earth map of {len(allr):,} georeferenced reference layers, "
            f"binned into 2-degree cells with circle area proportional to count. "
            f"Europe holds {eu_t:,} layers"
            + (f", up from {eu_b:,} before {LABEL.get(add_name, add_name)} was added." if add_name else "."))

    land = "\n".join(f'        <path d="{p}"/>' for p in land_paths())
    svg = f'''<svg viewBox="0 {vb_y:.0f} {W:.0f} {vb_h:.0f}" role="img" aria-label="{aria}">
      <g fill="none" stroke="currentColor" stroke-width=".5" opacity=".16">
{chr(10).join(f'        <path d="{p}"/>' for p in grat)}
      </g>
      <g fill="var(--sunken)" stroke="currentColor" stroke-width=".6" stroke-opacity=".38">
{land}
      </g>
      <g fill="var(--add)" fill-opacity=".5" stroke="var(--add)" stroke-width=".6" stroke-opacity=".85">
{circles(ca) if add_name else ""}
      </g>
      <g fill="var(--accent)" fill-opacity=".42" stroke="var(--accent)" stroke-width=".6" stroke-opacity=".8">
{circles(cb)}
      </g>
      <g font-family="IBM Plex Sans, sans-serif" font-size="11.5">
        <path d="{path(ex, ey, close=True)}" fill="none" stroke="currentColor" stroke-width="1.1" stroke-dasharray="4 3" opacity=".75"/>
        <text x="{ex.min():.0f}" y="{ey.min() - 9:.0f}" fill="currentColor" font-weight="600">{eu_label}</text>
        <text x="{ex.min():.0f}" y="{ey.max() + 16:.0f}" fill="var(--muted)">{eu_sub}</text>
      </g>
      <g font-family="IBM Plex Sans, sans-serif" font-size="11" fill="currentColor">{callouts}
      </g>
      <g font-family="IBM Plex Sans, sans-serif" font-size="11" fill="var(--muted)">
{size_svg}
        <text x="{sx:.0f}" y="{ly + 4:.0f}">layers per 2&#176; cell (area &#8733; count)</text>
      </g>
      <g font-family="IBM Plex Sans, sans-serif" font-size="11.5" fill="currentColor">
{key}
      </g>
    </svg>'''

    frag = os.path.join(HERE, "fig2_coverage_map.fragment.svg")
    open(frag, "w").write(svg + "\n")
    baked = svg
    for k, v in BAKE.items():
        baked = baked.replace(k, v)
    baked = baked.replace("<svg ", f'<svg xmlns="http://www.w3.org/2000/svg" '
                                   f'width="{W:.0f}" height="{vb_h:.0f}" ', 1)
    open(os.path.join(HERE, "fig2_coverage_map.svg"), "w").write(
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<!-- py_HB_texture coverage map. Land outline: Natural Earth 110 m, '
        'public domain. 2-degree cell counts only. -->\n' + baked + "\n")

    print(f"{'region':<15s} {base_name:>10s}" + (f" {'+' + add_name:>10s} {'added':>7s}" if add_name else ""))
    for k, box in REGIONS.items():
        b, t = in_region(base, box), in_region(allr, box)
        print(f"{k:<15s} {b:10,d}" + (f" {t:10,d} {t - b:+7,d}" if add_name else ""))
    print(f"{'georeferenced':<15s} {len(base):10,d}" + (f" {len(allr):10,d} {len(add):+7,d}" if add_name else ""))
    print(f"cells {len(cb)} base, {len(ca)} added; largest cell {nmax}")
    off = int((~allr.lat.between(LAT_S, LAT_N)).sum())
    print(f"off-map (outside {LAT_S:g}..{LAT_N:g} lat, not drawn): {off}")


if __name__ == "__main__":
    main()
