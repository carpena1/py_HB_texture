# Vertisols: converting retention curves with a shrinking bulk density

Swelling clays shrink as they dry, so their bulk density rises along the
retention curve. A curve made volumetric with one wet-state density
understates the water held at the dry end, and the tool then reads the soil as
too coarse. This folder holds the discussion of the correction and its effect,
applied to two vertisols: Laikipia (Kenya) and Tomejil (Andalusia, Spain).

| file | contents |
|---|---|
| `Vertisols_shrinkage_conversion_example.pdf` | the memo: method, a worked example for each soil, the comparison table, effect on texture and Ks, caveat |
| `vertisols_memo.html` | its source (print to PDF with headless Chrome) |
| `tomejil_shrinkage_vs_laikipia.png` / `.pdf` | the two soils side by side (curves, densities, how much the conversion changes θ) |
| `make_vertisol_figure.py` | draws that figure |
| `donana_shrinkage_check.png` / `.pdf`, `make_donana_check.py` | the contrast: the Doñana sands go through the same conversion unchanged |
| `test_caps.py`, `test_caps_results.csv` | texture and Ks of both soils as a new source under caps 1.4, 1.6, 1.75 and a per-sample cap; Laikipia's curves are the published ones, so its cap 1.6 row is the tool's own table |

The conversion: `w = θ / BD`, `BD(w) = min(cap, 1/(1/2.65 + w))`, `θ = w · BD(w)`
(normal shrinkage up to the shrinkage limit). The cap differs by soil:

- **Laikipia:** 1.6 g/cm³, as in the published dataset (Willard et al. 2026,
  HydroShare, doi:10.4211/hs.99080b48206e4b2c82e716516d463a61). The memo uses
  the published gravimetric curves and the core's KALRO density, and its
  converted curves equal the published volumetric series, which the tool's
  reference holds (`prepare_willard.py`).
- **Tomejil:** 1.75 g/cm³, the median oven-dry clod density of KSSL
  Vertisols with clay 45–65 % and 33 kPa density below 1.25, in
  `prepare_andalusia.py` (`shrinkage_bd`, `TOMEJIL_BD_MAX`).

Caps of 1.6 and 1.75 give the same texture results for both soils
(`test_caps.py`, results in `test_caps_results.csv`).

Both datasets are distributed (Laikipia published by the project team;
Andalusia with K. Vanderlinden's agreement), so this folder is public. The
scripts read data kept outside the repository: the Laikipia soils folder
downloaded from HydroShare into `data/willard_raw/hydroshare/`, and
`data/andalusia_raw/`. Run them from the repository root.
