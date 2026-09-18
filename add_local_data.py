"""Add your own lab-verified retention curves to the reference.

The largest gain the tool has shown comes from having a user's own data
source in the reference: +5 points of exact class, and Ks rank correlation
from 0.12 to 0.42. A sensor network that sends some of its sites to the lab
for texture (and ideally Ks) can feed those back here, so that later
predictions from the same network are matched against its own soils.

Input: one CSV, one row per retention point.
    sample_id, h_kpa, theta          required on every row
    sand, silt, clay                 required, % (per sample)
    ksat_cmh, depth_cm, bd,          optional (per sample); bd in g/cm3
    sample_type, site, lat, lon
Per-sample values may be repeated on every row of the sample or given once.
theta may be a fraction or % (detected per sample). sample_type defaults to
"undisturbed"; site groups samples from one location (defaults to sample_id).

Each sample is fitted with the tool's own van Genuchten routine and written to
data/local_reference.csv, which the default reference picks up automatically.
Adding a sample_id again replaces it. The file is git-ignored.

Usage:  python add_local_data.py mycurves.csv [--source NAME]
"""

import argparse
import os

import numpy as np
import pandas as pd

import swcc_texture as st

OUT = os.path.join(st.DATA_DIR, "local_reference.csv")
MIN_POINTS = 5          # four parameters need at least five points
MAX_RMSE = 0.03         # reject poor fits (m3/m3), as for the other tables
PER_SAMPLE = ("sand", "silt", "clay", "ksat_cmh", "depth_cm", "bd",
              "sample_type", "site", "lat", "lon")


def first(g, col):
    v = g[col].dropna() if col in g else pd.Series(dtype=object)
    return v.iloc[0] if len(v) else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("csv", help="retention points with per-sample properties")
    ap.add_argument("--source", default="local",
                    help="name of this data source, e.g. the sensor network")
    args = ap.parse_args()

    pts = pd.read_csv(args.csv, sep=None, engine="python")
    pts.columns = [c.strip().lower() for c in pts.columns]
    need = {"sample_id", "h_kpa", "theta", "sand", "silt", "clay"}
    if need - set(pts.columns):
        raise SystemExit(f"missing columns: {sorted(need - set(pts.columns))}")

    rows, rejected = [], []
    for sid, g in pts.groupby("sample_id", sort=False):
        g = g.dropna(subset=["h_kpa", "theta"])
        h, theta = g.h_kpa.to_numpy(float), g.theta.to_numpy(float)
        if np.any(theta > 1.0):
            theta = theta / 100.0
        p = {c: first(g, c) for c in PER_SAMPLE}
        frac = np.array([p["sand"], p["silt"], p["clay"]], float)
        if len(h) < MIN_POINTS:
            rejected.append((sid, f"{len(h)} points")); continue
        if not np.all(np.isfinite(frac)) or not 95 <= frac.sum() <= 105:
            rejected.append((sid, "sand + silt + clay not ~100 %")); continue
        sand, silt, clay = 100 * frac / frac.sum()
        try:
            popt, _ = st.fit_vg(h, theta)
        except Exception:
            rejected.append((sid, "fit failed")); continue
        tr, ts, la, ln1 = popt
        alpha, n = 10.0 ** la, 1.0 + 10.0 ** ln1
        rmse = float(np.sqrt(np.mean((st.vg_theta(h, tr, ts, alpha, n)
                                      - theta) ** 2)))
        if not np.isfinite(rmse) or rmse > MAX_RMSE:
            rejected.append((sid, f"fit RMSE {rmse:.3f}")); continue
        stype = str(p["sample_type"]).strip().lower()
        stype = stype if stype in st.SAMPLE_TYPE_CODE else "undisturbed"
        site = p["site"] if pd.notna(p["site"]) else sid
        rows.append(dict(
            layer_id=f"LOCAL_{args.source}_{sid}",
            profile_id=f"LOCAL_{args.source}_{site}",
            texture_class=st.usda_class(sand, silt, clay),
            alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
            sand=sand, silt=silt, clay=clay,
            ksat_cmh=p["ksat_cmh"] if pd.notna(p["ksat_cmh"]) and p["ksat_cmh"] > 0
            else np.nan,
            depth_cm=p["depth_cm"], bd=p["bd"], lat=p["lat"], lon=p["lon"],
            source_db=f"LOCAL_{args.source}", sample_type=stype,
            sample_type_source="supplied with the data", rmse=rmse,
            n_points=len(h)))

    new = pd.DataFrame(rows)
    old = pd.read_csv(OUT) if os.path.exists(OUT) else pd.DataFrame()
    replaced = int(old.layer_id.isin(new.get("layer_id", [])).sum()) if len(old) else 0
    if len(old) and len(new):
        old = old[~old.layer_id.isin(new.layer_id)]
    out = pd.concat([old, new], ignore_index=True)
    if len(out):
        out.to_csv(OUT, index=False)

    print(f"added {len(new) - replaced}, replaced {replaced}, rejected "
          f"{len(rejected)}; {OUT} now holds {len(out)} samples")
    for sid, why in rejected:
        print(f"  rejected {sid}: {why}")
    if len(new):
        print(f"  with Ks: {int(new.ksat_cmh.notna().sum())}  "
              f"classes: {new.texture_class.value_counts().to_dict()}")


if __name__ == "__main__":
    main()
