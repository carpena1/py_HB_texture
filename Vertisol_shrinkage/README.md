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
| `laikipia_shrinkage_example_WL_3LT07414.png` | the first, Laikipia-only example |
| `donana_shrinkage_check.png` / `.pdf`, `make_donana_check.py` | the contrast: the Doñana sands go through the same conversion unchanged |
| `Vertisol_calc.xlsx` | the memo's worked example for WL_3LT07414, step by step in a spreadsheet |
| `test_caps.py`, `test_caps_results.csv` | texture and Ks of both soils as a new source under caps 1.4, 1.6, 1.75 and a per-sample cap; Laikipia's curves come from `prepare_willard.scaled_curve()`, the conversion this memo develops |

The conversion: `w = θ / BD`, `BD(w) = min(cap, 1/(1/2.65 + w))`, `θ = w · BD(w)`
(normal shrinkage up to the shrinkage limit). The memo and figure use cap 1.75
g/cm³, the median oven-dry clod density of KSSL Vertisols with clay 45–65 % and
33 kPa density below 1.25; 1.6 (Tomejil's dry density) gives the same texture
results (`test_caps.py`, results in `test_caps_results.csv`).

What the tool's reference uses now:

- **Tomejil:** this conversion with cap 1.75, in `prepare_andalusia.py`
  (`shrinkage_bd`, `TOMEJIL_BD_MAX`).
- **Laikipia:** since 2026-10-16, the authors' own revised curves
  (L. Willard), converted with normal shrinkage up to 1.6 g/cm³, as the
  Scientific Data paper publishes them (`prepare_willard.py`, `REVISED_COL`).
  The conversion in this memo (cap 1.75 on the rescaled curve) is kept as
  `prepare_willard.scaled_curve()`, which `test_caps.py` uses. So the memo
  documents how the correction was developed, not the final Laikipia table.

Both datasets are now distributed (Laikipia, the project team's own data,
since 2026-10-31; Andalusia with K. Vanderlinden's agreement), so this folder
is public. The scripts read raw files that stay local (`data/willard_raw/`,
`data/andalusia_raw/`); run them from the repository root.
