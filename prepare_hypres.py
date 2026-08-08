"""Build a HYPRES reference table for use alongside GSHP.

HYPRES (Wosten et al. 1999) holds ~5,521 samples from 4,486 European soil
horizons. Per the JRC metadata it ships several tables, of which two matter
here:
    HYDRAULIC_PROPS  fitted Mualem-van Genuchten parameters per sample
    RAWRET           raw theta(h) data pairs (~197,000)
Either is enough: if fitted parameters are present they are used directly,
otherwise the raw pairs are fitted with the same routine the tool uses.

LICENCE -- READ FIRST
    HYPRES is distributed by ESDAC under terms that forbid passing the data to
    third parties. The output of this script is written to
    data/hypres_reference.csv, which is git-ignored, and MUST NOT be committed
    or published. Obtain the data under your own registration.

Because the exact column names of the distribution are not documented publicly,
this script matches columns by pattern and tells you what it found:

    python prepare_hypres.py --inspect <file-or-dir>     # show detected columns
    python prepare_hypres.py <file-or-dir>               # build the table

Accepts .csv, .txt, .xlsx and .mdb (the last needs `mdbtools`). Override any
detection with, e.g., --map alpha=ALFA --map n=ENN.
"""

import argparse
import glob
import io
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd

import swcc_texture as st

CM_TO_KPA = 0.0980665
OUT = "data/hypres_reference.csv"

# Field -> ordered list of regex patterns tried against lower-cased columns.
ALIASES = {
    "sample_id": [r"^sample(_?id)?$", r"^id$", r"^code$", r"^hor(izon)?_?id$"],
    "profile_id": [r"^profile(_?id)?$", r"^prof(_?id)?$", r"^site(_?id)?$",
                   r"^soil(_?id)?$"],
    "alpha": [r"^alp?h?a", r"^alfa", r"^vg_?a$"],
    "n": [r"^n$", r"^vg_?n$", r"^enn$"],
    "thetar": [r"^th?e?ta_?r$", r"^wcr$", r"^res.*water", r"^tr$"],
    "thetas": [r"^th?e?ta_?s$", r"^wcs$", r"^sat.*water", r"^ts$"],
    "sand": [r"^sand", r"^psand", r"^s_?perc"],
    "silt": [r"^silt", r"^psilt"],
    "clay": [r"^clay", r"^pclay"],
    "ksat": [r"^k_?sat", r"^ksat", r"^sat.*cond", r"^ks$"],
    "depth_top": [r"^(depth_?)?top", r"^hor.*top", r"^upper"],
    "depth_bot": [r"^(depth_?)?bot", r"^hor.*bot", r"^lower"],
    "depth": [r"^depth$", r"^mid_?depth"],
    "head": [r"^h$", r"^head", r"^pres", r"^suction", r"^tension", r"^pf$"],
    "theta": [r"^th?eta$", r"^wc$", r"^water_?content$", r"^vwc$"],
}


def read_any(path):
    """Read every table in a file or directory into {name: DataFrame}."""
    if os.path.isdir(path):
        out = {}
        for f in sorted(glob.glob(os.path.join(path, "*"))):
            out.update(read_any(f))
        return out
    ext = os.path.splitext(path)[1].lower()
    name = os.path.basename(path)
    try:
        if ext == ".mdb":
            tables = subprocess.run(["mdb-tables", "-1", path],
                                    capture_output=True, text=True,
                                    check=True).stdout.split()
            out = {}
            for t in tables:
                csv = subprocess.run(["mdb-export", path, t],
                                     capture_output=True, text=True).stdout
                if csv.strip():
                    out[f"{name}:{t}"] = pd.read_csv(io.StringIO(csv),
                                                     low_memory=False)
            return out
        if ext in (".xlsx", ".xls"):
            return {f"{name}:{s}": d for s, d
                    in pd.read_excel(path, sheet_name=None).items()}
        if ext in (".csv", ".txt", ".dat"):
            return {name: pd.read_csv(path, sep=None, engine="python",
                                      low_memory=False)}
    except Exception as e:                       # unreadable / not a table
        print(f"  ! could not read {name}: {e}", file=sys.stderr)
    return {}


def detect(df, overrides=None):
    """Map canonical field -> actual column name for one table."""
    found, used = {}, set()
    cols = {c: str(c).strip().lower() for c in df.columns}
    for field, pats in ALIASES.items():
        if overrides and field in overrides:
            if overrides[field] in df.columns:
                found[field] = overrides[field]
                used.add(overrides[field])
            continue
        for pat in pats:
            hit = [c for c, lc in cols.items()
                   if c not in used and re.match(pat, lc)]
            if hit:
                found[field] = hit[0]
                used.add(hit[0])
                break
    return found


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("source", help="HYPRES file or directory")
    ap.add_argument("--inspect", action="store_true",
                    help="only report the tables and columns detected")
    ap.add_argument("--map", action="append", default=[], metavar="FIELD=COL",
                    help="force a column mapping, e.g. --map alpha=ALFA")
    ap.add_argument("--alpha-units", choices=["1/cm", "1/kPa"], default="1/cm")
    ap.add_argument("--head-units", choices=["cm", "kPa", "pF"], default="cm")
    ap.add_argument("--ksat-units", choices=["cm/day", "cm/h"],
                    default="cm/day")
    args = ap.parse_args()

    overrides = dict(m.split("=", 1) for m in args.map)
    tables = read_any(args.source)
    if not tables:
        sys.exit(f"no readable tables found in {args.source}")

    print(f"found {len(tables)} table(s):")
    detected = {}
    for name, df in tables.items():
        f = detect(df, overrides)
        detected[name] = f
        print(f"  {name}  ({len(df)} rows, {len(df.columns)} cols)")
        print(f"      columns: {list(df.columns)[:14]}")
        if f:
            print(f"      detected: {f}")

    if args.inspect:
        print("\n(inspect only; re-run without --inspect to build, adding "
              "--map FIELD=COL for anything mis-detected)")
        return

    # A parameter table has vG parameters; a raw table has head/theta pairs.
    par = [(n, tables[n], d) for n, d in detected.items()
           if {"alpha", "n"} <= set(d)]
    raw = [(n, tables[n], d) for n, d in detected.items()
           if {"head", "theta"} <= set(d) and "alpha" not in d]
    tex = [(n, tables[n], d) for n, d in detected.items()
           if {"sand", "silt", "clay"} <= set(d)]

    if par:
        name, df, d = par[0]
        print(f"\nusing fitted vG parameters from '{name}'")
        rec = pd.DataFrame({
            "alpha_kpa": pd.to_numeric(df[d["alpha"]], errors="coerce"),
            "n": pd.to_numeric(df[d["n"]], errors="coerce"),
            "thetar": pd.to_numeric(df[d["thetar"]], errors="coerce")
            if "thetar" in d else np.nan,
            "thetas": pd.to_numeric(df[d["thetas"]], errors="coerce")
            if "thetas" in d else np.nan,
        })
        if args.alpha_units == "1/cm":
            rec["alpha_kpa"] /= CM_TO_KPA
        key = d.get("sample_id")
        rec["layer_id"] = ("HYPRES_" + df[key].astype(str)) if key else \
            [f"HYPRES_{i}" for i in range(len(df))]
        rec["profile_id"] = ("HYPRES_P" + df[d["profile_id"]].astype(str)) \
            if "profile_id" in d else rec["layer_id"]
        for f, col in (("ksat_cmh", "ksat"), ("depth_cm", "depth")):
            rec[f] = pd.to_numeric(df[d[col]], errors="coerce") \
                if col in d else np.nan
        if "ksat" in d and args.ksat_units == "cm/day":
            rec["ksat_cmh"] /= 24.0
        if "depth" not in d and {"depth_top", "depth_bot"} <= set(d):
            rec["depth_cm"] = (pd.to_numeric(df[d["depth_top"]],
                                             errors="coerce") +
                               pd.to_numeric(df[d["depth_bot"]],
                                             errors="coerce")) / 2.0
        src = df
    elif raw:
        name, df, d = raw[0]
        print(f"\nfitting vG to raw theta(h) pairs from '{name}'")
        key = d.get("sample_id") or df.columns[0]
        h = pd.to_numeric(df[d["head"]], errors="coerce")
        if args.head_units == "pF":
            h = 10.0 ** h
        h = h * (CM_TO_KPA if args.head_units in ("cm", "pF") else 1.0)
        work = pd.DataFrame({"id": df[key], "h": h,
                             "theta": pd.to_numeric(df[d["theta"]],
                                                    errors="coerce")}).dropna()
        rows = []
        for sid, g in work.groupby("id"):
            if len(g) < 5:
                continue
            try:
                popt, _ = st.fit_vg(g.h.to_numpy(), g.theta.to_numpy())
            except Exception:
                continue
            tr, ts, la, ln1 = popt
            rows.append(dict(layer_id=f"HYPRES_{sid}",
                             profile_id=f"HYPRES_{sid}", thetar=tr, thetas=ts,
                             alpha_kpa=10.0 ** la, n=1.0 + 10.0 ** ln1,
                             ksat_cmh=np.nan, depth_cm=np.nan))
        rec = pd.DataFrame(rows)
        src = None
    else:
        sys.exit("no table with van Genuchten parameters or theta(h) pairs "
                 "was recognised; re-run with --inspect and use --map")

    # Attach particle fractions, from the same or another table.
    if tex:
        tname, tdf, td = tex[0]
        key = td.get("sample_id")
        f = pd.DataFrame({
            "sand": pd.to_numeric(tdf[td["sand"]], errors="coerce"),
            "silt": pd.to_numeric(tdf[td["silt"]], errors="coerce"),
            "clay": pd.to_numeric(tdf[td["clay"]], errors="coerce")})
        f["layer_id"] = ("HYPRES_" + tdf[key].astype(str)) if key else \
            [f"HYPRES_{i}" for i in range(len(tdf))]
        rec = rec.merge(f, on="layer_id", how="left")
        print(f"attached particle fractions from '{tname}'")
    else:
        sys.exit("no sand/silt/clay columns found; texture is required")

    tot = rec[["sand", "silt", "clay"]].sum(axis=1)
    ok = tot.between(95, 105)
    for c in ("sand", "silt", "clay"):
        rec[c] = (rec[c] / tot * 100).where(ok)
    rec["texture_class"] = [
        st.usda_class(s, si, c) if np.isfinite(s) else None
        for s, si, c in zip(rec["sand"], rec["silt"], rec["clay"])]

    rec = rec[rec.texture_class.notna() & (rec.n > 1.0) & (rec.alpha_kpa > 0)]
    rec = rec[["layer_id", "profile_id", "texture_class", "alpha_kpa", "n",
               "thetar", "thetas", "sand", "silt", "clay", "ksat_cmh",
               "depth_cm"]]
    rec.to_csv(OUT, index=False)
    print(f"\nwrote {OUT}: {len(rec)} soils")
    print(rec.texture_class.value_counts())

    gshp = pd.read_csv(st.REFERENCE_CSV)
    print("\nsanity check -- class-median alpha [kPa^-1] / n, HYPRES vs GSHP")
    print("(a ~10x offset in alpha means --alpha-units is wrong)")
    for cls in ["sand", "sandy loam", "silt loam", "clay"]:
        a, b = rec[rec.texture_class == cls], gshp[gshp.texture_class == cls]
        if len(a) and len(b):
            print(f"  {cls:<12s} {a.alpha_kpa.median():7.3f} vs "
                  f"{b.alpha_kpa.median():7.3f}   n {a.n.median():5.2f} vs "
                  f"{b.n.median():5.2f}   ({len(a)} vs {len(b)} soils)")
    print("\nREMINDER: data/hypres_reference.csv is git-ignored. Do not commit "
          "or redistribute it.")


if __name__ == "__main__":
    main()
