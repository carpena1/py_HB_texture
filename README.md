# py_HB_texture — soil texture & Ks from the water characteristic curve

Estimates USDA texture class, particle fractions (% sand/silt/clay) and
saturated hydraulic conductivity Ks (cm/h), with uncertainty ranges, from
measured SWCC data (h in kPa, theta as fraction or % — auto-detected).

Method (see Problem_description.docx): fit van Genuchten parameters
(Mualem m = 1-1/n), then distance-weighted k-nearest-neighbor inference in
the GSHP database (Gupta et al. 2022, doi:10.5281/zenodo.6640246;
9,996 quality-filtered layers). Votes use inverse class-frequency weighting
(uniform prior over the 12 USDA classes). Fit uncertainty is propagated by
Monte Carlo sampling of the fit covariance.

## Setup

    python3 -m venv .venv
    .venv/bin/pip install -r requirements.txt

Calling `.venv/bin/python` directly (as in Usage below) runs inside the
virtual environment without activating it, so no extra step is needed. If you
prefer to type `python` instead of the full path, activate the environment
first with `source .venv/bin/activate` (`deactivate` to exit).

The derived reference table `data/gshp_reference.csv` is included, so the tool
runs out of the box. The 80 MB raw GSHP file is not committed; to regenerate
the reference, download it from Zenodo (record 6640246) into `data/` and run
`prepare_gshp.py`.

## Usage

    .venv/bin/python swcc_texture.py mydata.csv [--json out.json]

`mydata.csv`: two columns (comma or whitespace), h_kPa and theta.
Or from Python: `swcc_texture.estimate(h, theta)`.

## Files

- `swcc_texture.py` — the tool (fit + inference + CLI)
- `prepare_gshp.py` — one-time: raw GSHP csv -> `data/gshp_reference.csv`
- `verify_carsel_parrish.py` — independent verification vs Carsel & Parrish (1988)
- `verify_rosetta.py` — independent verification vs ROSETTA class means
- `verify_gshp.py` — leave-one-out verification vs source GSHP soils
- `verify_groups.py` — accuracy at the 4 aggregated texture groups
- `prepare_unsoda.py` — one-time: UNSODA 2.0 .mdb -> `data/unsoda_reference.csv`
- `verify_variants.py` — reference/covariate variants on the 4 benchmarks
- `verify_variants_large.py` — same variants, large paired leave-one-out test
- `data/gshp_reference.csv` — the distributed reduced GSHP reference table
  (~10 k curated layers), derived from the raw GSHP file via `prepare_gshp.py`.
  The 80 MB raw file (`WRC_dataset_surya_et_al_2021_final.csv`, Zenodo 6640246)
  is not included; download it only if you want to rebuild the reference.
- `data/unsoda_reference.csv` — 588 UNSODA 2.0 soils in the same schema, used
  by the merged-reference variant (raw .mdb not committed; rebuild with
  `prepare_unsoda.py`, which needs `mdbtools`)
- `Testing/` — per-class test curves: `<i>_<code>_CnP.csv` (Carsel & Parrish
  1988) and `<i>_<code>_GHSP.csv` (a representative GSHP soil for each class)

## Verification vs. Carsel & Parrish (1988) — corrected loam alpha (2026-07-21)

Synthetic curves from the Carsel & Parrish class-typical van Genuchten
parameters (the per-class curves are shipped as `Testing/<i>_<code>_CnP.csv`).
Texture codes: S sand, LS loamy sand, SL sandy loam, L loam, Si silt,
SiL silt loam, SCL sandy clay loam, CL clay loam, SiCL silty clay loam,
SC sandy clay, SiC silty clay, C clay.

Per-class accuracy table (`verify_carsel_parrish.py`). `n_ref` = curated GSHP
layers backing texture/fractions for that class; `n_ks` = subset with a measured
Ksat; `nb` = mean number of the 30 nearest neighbors carrying a measured Ksat
(the support behind each Ks estimate).

| True class | n_ref | n_ks | Predicted | P(true) | Rank | Ks pred | Ks true | nb | Ks in 5-95% |
|---|---|---|---|---|---|---|---|---|---|
| sand (S) | 3930 | 3735 | **sand (S)** | 1.00 | 1 | 4.55 | 29.7 | 26 | no |
| loamy sand (LS) | 805 | 581 | sand (S) | 0.09 | 2 | 3.62 | 14.6 | 24 | no |
| sandy loam (SL) | 1437 | 785 | **sandy loam (SL)** | 0.39 | 1 | 1.86 | 4.42 | 25 | yes |
| loam (L) | 629 | 103 | loamy sand (LS) | 0.00 | 12 | 0.33 | 1.04 | 23 | yes |
| silt (Si) | 32 | 10 | silt loam (SiL) | 0.00 | 12 | 0.059 | 0.25 | 16 | yes |
| silt loam (SiL) | 665 | 191 | **silt loam (SiL)** | 0.55 | 1 | 0.14 | 0.45 | 15 | yes |
| sandy clay loam (SCL) | 918 | 660 | **sandy clay loam (SCL)** | 0.50 | 1 | 0.20 | 1.31 | 14 | yes |
| clay loam (CL) | 291 | 60 | silty clay (SiC) | 0.09 | 6 | 0.053 | 0.26 | 21 | yes |
| silty clay loam (SiCL) | 284 | 49 | silty clay (SiC) | 0.14 | 3 | 0.18 | 0.07 | 11 | yes |
| sandy clay (SC) | 161 | 108 | sandy loam (SL) | 0.00 | 12 | 0.023 | 0.12 | 25 | yes |
| silty clay (SiC) | 131 | 39 | sandy clay loam (SCL) | 0.14 | 3 | 0.081 | 0.02 | 25 | yes |
| clay (C) | 713 | 323 | sandy clay loam (SCL) | 0.02 | 8 | 0.10 | 0.20 | 23 | yes |
| **Total** | **9996** | **6644** | | | | | | | |

- vG fit stage: exact recovery of all 12 generating parameter sets.
- Texture class: 4/12 exact, 5/12 top-2. Development-set (GSHP leave-one-out)
  ceiling is ~41 % macro recall even using true vG parameters — the
  vG -> texture mapping is inherently ambiguous.
- Fractions: mean |error| vs USDA class centroid 14.3 %.
- Ks: benchmark value inside the 5-95 % range for 10/12 classes;
  typical factor ~3.8. GSHP lab sands are much slower than the idealized
  Carsel & Parrish sand (4.6 vs 29.7 cm/h).
- Curated reference: 9,996 layers back the texture/fraction inference; only
  6,644 carry a measured Ksat, so Ks rests on a smaller, sand-skewed subset
  (per class: loam 103, clay loam 60, silty clay loam 49 with Ksat). The
  verification table reports per-class n_ref/n_ks and 'nb', the mean number of
  the 30 nearest neighbors that had a measured Ksat behind each Ks estimate;
  the CLI prints the same neighbor count per run.
- Dominant error source: the benchmark is built from Carsel & Parrish (1988)
  class-typical parameters, which differ systematically from GSHP's empirical
  class-conditional parameters (e.g. GSHP sand alpha median 0.25 vs Carsel
  1.45 kPa^-1; GSHP clay n median 1.36, thetar 0.22 vs Carsel 1.09, 0.068).
  The tool reproduces GSHP's mapping; where Carsel and GSHP disagree, the
  benchmark penalizes the tool regardless of method.

## Verification vs ROSETTA class means (Schaap et al. 2001; verify_rosetta.py)

Second benchmark: synthetic curves from the ROSETTA H1 class-average vG
parameters (a more modern, data-derived parameterization than Carsel).
Predicted fractions are compared to each class's USDA centroid (as for
Carsel); `nb` is the mean number of the 30 neighbors with a measured Ksat.

| True class | Predicted | P(true) | Rank | Ks pred | Ks true | nb | Ks in 5-95% |
|---|---|---|---|---|---|---|---|
| sand (S) | **sand (S)** | 1.00 | 1 | 1.86 | 26.8 | 30 | no |
| loamy sand (LS) | sand (S) | 0.32 | 2 | 0.59 | 4.38 | 24 | yes |
| sandy loam (SL) | loamy sand (LS) | 0.11 | 4 | 0.24 | 1.60 | 26 | no |
| loam (L) | silt loam (SiL) | 0.13 | 2 | 1.73 | 0.50 | 18 | yes |
| silt (Si) | clay loam (CL) | 0.00 | 12 | 0.41 | 1.82 | 2 | no |
| silt loam (SiL) | **silt loam (SiL)** | 0.28 | 1 | 0.036 | 0.76 | 1 | no |
| sandy clay loam (SCL) | sandy loam (SL) | 0.10 | 6 | 0.079 | 0.55 | 24 | yes |
| clay loam (CL) | silt loam (SiL) | 0.04 | 7 | 0.45 | 0.34 | 15 | yes |
| silty clay loam (SiCL) | clay loam (CL) | 0.10 | 3 | 2.42 | 0.46 | 1 | no |
| sandy clay (SC) | loam (L) | 0.00 | 12 | 2.09 | 0.47 | 22 | yes |
| silty clay (SiC) | silty clay loam (SiCL) | 0.00 | 12 | 0.045 | 0.40 | 11 | yes |
| clay (C) | silt loam (SiL) | 0.00 | 12 | 0.042 | 0.62 | 17 | yes |

- vG fit stage: exact recovery of all 12 ROSETTA parameter sets (RMSE 0).
- Texture class: 2/12 exact, 4/12 top-2.
- Fractions: mean |error| vs centroid 12.9 % (centroid in 5-95 % range 7/12).
- Ks: true value in 5-95 % range 7/12; typical factor ~6.6.
- This scores *below* the Carsel benchmark, and the reason is instructive:
  the three databases disagree on where each class sits in vG space. Class-mean
  alpha [kPa^-1] for clay is 0.08 (Carsel) / 0.15 (ROSETTA) / 0.41 (GSHP) — a
  5x spread. Comparing ROSETTA class-means directly to GSHP class-medians
  (no kNN, cleanest possible test) still recovers only 4/12 classes, so the
  low score reflects database disagreement, not a tool defect. ROSETTA's
  classes are also more compressed in alpha (0.05-0.36 vs Carsel 0.05-1.45),
  so its fine classes collapse toward the GSHP loam/clay-loam center.

Realistic accuracy is bounded by GSHP's own internal consistency (~41 %
leave-one-out macro recall); no lookup tool exceeds that against a foreign
benchmark. The tool therefore reports full class probabilities and Ks/fraction
ranges rather than a single hard answer.

## Verification vs source GSHP, 12 representative soils (verify_gshp.py)

Final, in-distribution check. For each USDA class we take the real GSHP soil
closest to that class's median in vG space (the soils shipped as
`Testing/<i>_<code>_GHSP.csv`), then classify it with that exact soil **removed
from the reference** (leave-one-out), so it cannot match itself. Because these
are real soils, predictions are checked against their **measured** sand/silt/clay
and Ks — not class centroids.

| True class | Predicted | P(true) | Rank | s/si/cl pred | s/si/cl true | Ks pred | Ks true | nb | Ks in 5-95% |
|---|---|---|---|---|---|---|---|---|---|
| sand (S) | **sand (S)** | 1.00 | 1 | 94/4/2 | 96/3/1 | 0.85 | 0.69 | 30 | yes |
| loamy sand (LS) | **loamy sand (LS)** | 0.35 | 1 | 83/7/9 | 86/7/7 | 0.26 | 0.25 | 22 | yes |
| sandy loam (SL) | loamy sand (LS) | 0.28 | 3 | 81/12/7 | 61/31/8 | 0.099 | 0.058 | 23 | yes |
| loam (L) | silt loam (SiL) | 0.18 | 4 | 54/32/15 | 50/33/17 | 0.092 | 0.048 | 22 | yes |
| silt (Si) | **silt (Si)** | 0.93 | 1 | 4/85/10 | 6/84/10 | 0.061 | 0.37 | 22 | yes |
| silt loam (SiL) | silt (Si) | 0.02 | 7 | 21/53/25 | 26/64/10 | 0.036 | 2.25 | 2 | no |
| sandy clay loam (SCL) | **sandy clay loam (SCL)** | 0.45 | 1 | 67/14/19 | 71/6/23 | 0.035 | 0.039 | 27 | yes |
| clay loam (CL) | **clay loam (CL)** | 0.32 | 1 | 18/48/34 | 23/44/33 | 7.84 | 70.5 | 7 | no |
| silty clay loam (SiCL) | **silty clay loam (SiCL)** | 0.41 | 1 | 22/48/30 | 7/54/39 | 3.79 | 4.54 | 1 | no |
| sandy clay (SC) | sandy clay loam (SCL) | 0.00 | 12 | 55/22/23 | 59/4/36 | 0.032 | 0.036 | 15 | yes |
| silty clay (SiC) | silty clay loam (SiCL) | 0.06 | 5 | 28/36/36 | 11/46/43 | 0.095 | 1.11 | 4 | no |
| clay (C) | **clay (C)** | 0.40 | 1 | 30/21/49 | 45/9/46 | 4.75 | 0.43 | 11 | no |

- Texture class: **7/12 exact, 7/12 top-2**.
- Fractions: mean |error| vs the soils' **measured** fractions **6.8 %**;
  measured value inside the 5-95 % range for 11/12 soils.
- Ks: measured value inside the 5-95 % range for 7/12; typical factor ~3.4.

This is the fair "does it work on real soils it wasn't allowed to memorize"
test, and it is the strongest of the three (7/12 vs 4/12 Carsel, 2/12 ROSETTA)
precisely because the inputs come from the same population as the reference.

## Accuracy at aggregated texture groups (verify_groups.py)

The 12 USDA classes collapsed into four broad groups; a prediction is scored
correct if it lands in the true class's group (a same-group miss is a
"near miss"). This re-scores the three verifications above — the tool itself is
unchanged; only the evaluation is coarsened.

| Group | Classes |
|---|---|
| Sandy | sand, loamy sand |
| Loamy | sandy loam, loam, sandy clay loam, clay loam |
| Silty | silt, silt loam, silty clay loam |
| Clayey | sandy clay, silty clay, clay |

| Benchmark | Class-level (exact) | Group-level |
|---|---|---|
| Carsel & Parrish (1988) | 4/12 | 6/12 (50 %) |
| ROSETTA class means | 2/12 | 4/12 (33 %) |
| GSHP leave-one-out | 7/12 | 8/12 (67 %) |

Grouping helps least where it might be expected to help most: the **Clayey**
group is the persistent failure across all three benchmarks — real/typical clays
(sandy clay, silty clay, clay) are pulled into the Loamy or Silty clay-loam
family, because in van Genuchten space the fine classes overlap heavily and GSHP
has few of them. The other recurring near-misses are the adjacent
sandy-loam/loamy-sand and loam/silt-loam boundaries. On the in-distribution GSHP
test the four group misses are all single-step, adjacent-group confusions.

## References

- Carsel, R.F. and Parrish, R.S. (1988). Developing joint probability
  distributions of soil water retention characteristics. *Water Resources
  Research* 24(5):755–769. doi:10.1029/WR024i005p00755.
  (Source of the class-typical van Genuchten parameters used in the
  `Testing/<i>_<code>_CnP.csv` benchmark curves.)
- Gupta, S., Papritz, A., Lehmann, P., Hengl, T., Bonetti, S. and Or, D.
  (2022). Global Soil Hydraulic Properties dataset based on legacy site
  observations and robust parameterization. *Scientific Data* 9:444.
  doi:10.1038/s41597-022-01481-5. Data: doi:10.5281/zenodo.6640246 (CC BY 4.0).
  (The GSHP reference database used for inference.)
- Schaap, M.G., Leij, F.J. and van Genuchten, M.Th. (2001). ROSETTA: a computer
  program for estimating soil hydraulic parameters with hierarchical
  pedotransfer functions. *Journal of Hydrology* 251:163–176.
  (Source of the ROSETTA class-mean parameters used in `verify_rosetta.py`.)
- van Genuchten, M.Th. (1980). A closed-form equation for predicting the
  hydraulic conductivity of unsaturated soils. *Soil Science Society of
  America Journal* 44:892–898.

## License

Released under the Creative Commons Attribution 4.0 International license
(CC BY 4.0); full text in `LICENSE`. You may share and adapt this work,
including commercially, provided you give appropriate credit. Suggested
attribution: "R. Muñoz-Carpena, py_HB_texture, https://github.com/carpena1/py_HB_texture".

The bundled GSHP reference data (`data/gshp_reference.csv`, derived from
Gupta et al. 2022) is itself CC BY 4.0 — cite Gupta et al. and Zenodo record
6640246 when reusing it.
