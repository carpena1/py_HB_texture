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

Optional flags (see "Improving accuracy" and "Method levers" below):

    --model hybrid             predict the class with a gradient-boosted
                               classifier instead of neighbour voting:
                               +3.2 pp (p=0.005), +5.4 pp against an unseen
                               laboratory. RECOMMENDED. Needs scikit-learn.
    --tau 0.75                 temper the class prior; better macro-F1 and
                               far better precision for rare classes
    --depth 25                 sample mid-depth in cm; +3.5 pp on average, but
                               HARMFUL for European soils (see below)
    --reference merged         add the 588 UNSODA 2.0 soils to the reference

## Files

- `swcc_texture.py` — the tool (fit + inference + CLI)
- `prepare_gshp.py` — one-time: raw GSHP csv -> `data/gshp_reference.csv`
- `verify_carsel_parrish.py` — independent verification vs Carsel & Parrish (1988)
- `verify_rosetta.py` — independent verification vs ROSETTA class means
- `verify_gshp.py` — leave-one-out verification vs source GSHP soils
- `verify_groups.py` — accuracy at the 4 aggregated texture groups
- `prepare_unsoda.py` — one-time: UNSODA 2.0 .mdb -> `data/unsoda_reference.csv`
- `prepare_sdb.py` — one-time: sDB .mat -> `data/sdb_reference.csv` (CC-BY-3.0)
- `verify_variants.py` — reference/covariate variants on the 4 benchmarks
- `verify_variants_large.py` — same variants, large paired leave-one-out test
- `verify_region.py` — how much region-matched reference data is worth
- `verify_hypres.py` — European benchmark from the HYPRES class PTFs
- `prepare_hypres.py` — loader kept in case sample-level HYPRES is ever released
- `ksat_metrics.py` — Ksat scoring shared by every verification script
- `verify_by_class.py` — per-class/per-group accuracy, with and without depth
- `verify_skill.py` — is each class acceptable in absolute terms?
- `prepare_kssl.py` — one-time: NCSS/KSSL SQLite -> `data/kssl_reference.csv`
- `verify_kssl.py` — does merging KSSL help?
- `verify_source_blocked.py` — accuracy against an unseen laboratory
- `verify_tau.py` — tuning the class-prior exponent
- `verify_features.py` — vG parameters vs curve-space matching
- `verify_gbm.py` — gradient boosting vs kNN on matched folds
- `verify_ceiling.py` — model-free Cover & Hart bound on achievable accuracy
- `verify_hybrid.py` — end-to-end check of the shipped hybrid model
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

## Improving accuracy: reference database and covariates

Two routes to better classification were tested (`verify_variants.py`,
`verify_variants_large.py`):

1. **Merge a second database.** HYPRES was the original candidate. The database
   as *built* holds exactly what this tool needs (the JRC metadata documents
   RAWRET/RAWK tables with ~197,000 theta(h) pairs and per-sample Mualem-van
   Genuchten parameters for 5,521 samples from 4,486 horizons), but the
   released package contains none of it — see "Outcome: HYPRES cannot be
   merged" below. **UNSODA 2.0** (Nemes et al. 2001) was used instead: a
   genuinely redistributable, sample-level database. 588 of its 790 soils
   yielded usable retention curves, fitted here with the same vG routine
   (median RMSE 0.0065).
   Note that GSHP already contains 218 UNSODA-sourced layers, so the merge is
   partly redundant — a likely reason its benefit is small.
2. **Add a covariate.** Sample **mid-depth** (98 % coverage in GSHP) is used as
   a fifth, equally weighted kNN feature. *Organic matter was dropped*: it is
   measured for only 15 % of curated GSHP layers, so using it would shrink the
   reference to a small, sand-skewed subset.

Because the Carsel and ROSETTA benchmarks are class-typical parameter sets with
no depth attached, covariates can only be scored on real soils. The decisive
test is therefore a **large paired leave-one-out** over 2,680 real soils
(2,121 GSHP + 559 UNSODA, stratified by class), each predicted with itself
removed from the reference and compared variant-vs-baseline on the *same* soils
(McNemar).

93 % of GSHP layers belong to multi-layer profiles (median 4 horizons), and
sibling horizons share both depth and texture — so a naive test could let depth
find siblings instead of genuinely similar soils. The run is therefore repeated
dropping the target's **whole profile**. The depth gain survives intact
(+3.5 pp vs +3.2 pp), confirming real pedological signal rather than site
leakage.

Profile-level leave-one-out, n = 2,680 paired soils:

| Variant | Class | vs base | p | Group | vs base | p |
|---|---|---|---|---|---|---|
| GSHP (baseline) | 39.8 % | — | — | 62.3 % | — | — |
| GSHP + depth | 43.3 % | **+3.5 pp** | <0.001 | 62.8 % | +0.6 pp | 0.52 |
| GSHP + UNSODA | 40.5 % | +0.7 pp | 0.052 | 62.7 % | +0.4 pp | 0.29 |
| GSHP + UNSODA + depth | 44.4 % | **+4.7 pp** | <0.001 | 63.6 % | +1.3 pp | 0.12 |

Conclusions:

- **Depth is worth using — but only where the reference is regionally dense.**
  Globally it gives +3.5 pp exact-class accuracy, highly significant and robust
  to profile-level leave-one-out. That global average, however, hides a
  regional reversal; see "Depth has opposite signs in dense and sparse regions"
  below before using `--depth`.
- **Merging UNSODA gives little**: +0.7 pp, at the edge of significance
  (p = 0.052). It adds sample diversity at no cost, so it is available via
  `--reference merged`, but it is not the default and should not be expected to
  change results much.
- **Both gains are at the class level only.** Group-level accuracy barely moves
  (+0.6 pp, p = 0.52), i.e. depth sharpens distinctions *within* a texture group
  (sand vs loamy sand) but does not fix the coarse Clayey-group confusions that
  dominate the group-level error.
- The 12-soil benchmark tables are too small to see these effects: there, depth
  appears to *hurt* (16/24 -> 15/24). A one-soil difference at n = 12 is noise;
  the n = 2,680 paired test is the one to trust.

### Depth has opposite signs in dense and sparse regions

The global +3.5 pp for `--depth` is an average that conceals a reversal. Tested
on European targets against a **class-matched** non-European control of equal
size (653 each, same reference, profile-level leave-one-out):

| Targets | No depth | With depth | Effect |
|---|---|---|---|
| European | 34.8 % | 27.7 % | **-7.0 pp** (p < 0.001) |
| Non-European (class-matched) | 35.8 % | 41.5 % | **+5.7 pp** (p = 0.003) |

Same reference, same class mix, opposite sign. Matching the class distribution
rules out texture composition as the explanation.

The likely mechanism is sparsity. Depth adds a fifth dimension to the matching
space, and Europe contributes only 778 of 9,996 layers. Where the reference is
regionally dense, depth genuinely narrows the neighbourhood onto pedologically
similar soils; where it is thin, the extra dimension pulls in soils that match
on depth but come from an entirely different setting, displacing the few
regionally appropriate neighbours. This is the curse of dimensionality biting
exactly where coverage is weakest.

Two consequences:

- **Do not use `--depth` for European soils** — it costs about 7 points there.
  Use it where the reference is dense for the region in question.
- It sharpens the case for European reference data: the coverage gap does not
  merely cost accuracy, it inverts the value of an otherwise useful covariate.
  On top of `--depth`, adding sDB's European horizons recovers +1.7 pp
  (p = 0.019), whereas without depth the same horizons do nothing (+0.3 pp,
  p = 0.83).

Caveat: "European" is a bounding box, so the effect is regional in the sense of
being associated with that subset of the database; it may be partly confounded
with the particular source datasets and measurement methods that dominate
there. The sparsity explanation is a hypothesis consistent with the numbers,
not something these runs prove directly.

### Would HYPRES help? Evidence from a regional hold-out (verify_region.py)

HYPRES cannot be tested without the data, but its *mechanism* can: it would add
~4,486 European horizons to a reference holding only 778 (7.8 % of curated
GSHP). So we measured what region-matched reference data is worth, by removing
European soils from the reference and re-classifying European targets. Three
conditions, target's own profile always excluded, n = 778 European soils:

| Condition | Reference | Class | Group |
|---|---|---|---|
| A full | 9,996 | 36.2 % | 60.2 % |
| B Europe removed | 9,218 | **27.5 %** | **52.3 %** |
| C 778 random non-European removed | 9,218 | 36.8 % | 61.7 % |

B vs C is the informative contrast — same reference size, only the region match
differs: **-9.3 pp class and -9.4 pp group (both p < 0.001)**. Simply having a
smaller reference costs nothing (A vs C: +0.5 pp, p = 0.50). A non-European
control group is unaffected by removing European data (+0.0 pp, p = 1.00),
confirming the effect is regional and not an artifact of reference size.

Two further points make the case: European soils are currently classified far
worse than non-European ones (36.2 % vs 52.8 % exact class — a 16.6 pp deficit
consistent with under-representation), and HYPRES would raise European coverage
roughly six-fold. This is a much larger prospective gain than either the depth
covariate (+3.5 pp) or the UNSODA merge (+0.7 pp) delivered — though it is an
estimate of the mechanism, not a measurement of HYPRES itself.

### Outcome: HYPRES cannot be merged (data was never released)

The ESDB v2.0 distribution was obtained (ESDAC access ID 132133, August 2026).
`Hypres.zip` contains three documents and no data tables. `HYPRES_Readme.doc`
states that the release contains the metadata and "the actual functions and
function Parameters", and that it

> does not contain ... HYPRES project source data and results as no agreement
> has been reached with the participating institutions regarding their
> distribution.

So the sample-level tables (RAWRET / RAWK / HYDRAULIC_PROPS) were never
distributed. This is not an ESDAC licensing decision that could be waived on
request: rights rest with the 20 contributing institutions from 12 countries,
and no distribution agreement was reached. The JRC metadata describes the
*internal* structure of the database as built, which is why it lists those
tables; the released package is a different, much smaller thing.

The companion Soil Profile Analytical Database (SPADBE) was also checked, since
it carries water-retention fields at -1, -10, -100 and -1500 kPa. Of its
measured-profile horizons, exactly **2** have both a texture analysis and >= 4
real retention values (the rest are the -999 missing-data sentinel), so it
cannot contribute to the reference either.

What HYPRES *does* provide is 11 class-average Mualem-van Genuchten parameter
sets (5 topsoil, 5 subsoil, 1 organic) plus continuous PTFs. Eleven points
among ~10,000 cannot move a kNN reference, but they make a useful **European
benchmark** — see below. The continuous PTFs run texture -> vG, the opposite of
this tool's direction, so they are not usable here either.

### Verification vs HYPRES class PTFs (verify_hypres.py)

A European counterpart to the US-centric Carsel and ROSETTA benchmarks.
Synthetic curves are generated from the 10 mineral class PTFs; since HYPRES
uses FAO/SGDBE classes rather than USDA ones, the tool's predicted
sand/silt/clay is scored against the FAO class envelope (coarse: clay < 18 and
sand > 65; medium; medium fine: clay < 35 and sand < 15; fine: 35 <= clay < 60;
very fine: clay >= 60).

- FAO texture class correct: **5/10**, and **6/10** when `--depth` is supplied
  (topsoil 15 cm, subsoil 60 cm) — independent corroboration that the depth
  covariate helps.
- Ks: HYPRES value inside the 5-95 % range for 9/10 classes, but the median is
  off by a typical factor of 10 — worse than on the other benchmarks.
- The failures are concentrated at the fine end: the "very fine" class
  (clay >= 60 %) is never recovered, with the tool predicting 41-48 % clay. This
  is the same fine-texture weakness the Carsel, ROSETTA and group-level results
  all show, now confirmed on European data.

### Survey: what European hydraulic data is actually obtainable

Prompted by the ~9 pp regional deficit, the candidate sources were surveyed.
The requirement is strict: *sample-level* data pairing a retention curve (or
fitted vG parameters) **with** a particle-size analysis. That combination is
rare.

| Source | Size | Status |
|---|---|---|
| GSHP | 9,996 curated (778 European) | **in use**, CC BY 4.0 |
| UNSODA 2.0 | 588 usable | **in use**, public (Ag Data Commons) |
| HYPRES | 4,486 horizons | **unavailable** — raw data never released |
| EU-HYDI | >18,000 samples, 29 institutions, 18 countries | **unavailable** — see below |
| sDB (Vereecken et al. 2017) | 182 horizons, 165 usable | **open (CC-BY-3.0)**, included — but ~88 % already in GSHP, so no measurable gain |
| SoilKsatDB | 13,258 Ksat + 11,584 texture | open (Zenodo) but **no retention curves**, so it cannot join a vG-matched reference; its Ksat also overlaps GSHP (same authors) |
| HYDROS | 173 samples | too small to matter |
| HYBRAS | 445 sites | Brazilian; already contributes 814 layers to GSHP |

EU-HYDI is the only remaining large European candidate, and its report
(EUR 26053 EN, section 2 "Notice on data access") is explicit:

> All raw data contained in the EU-HYDI database are accessible only to the
> contributing participants and the EU-HYDI coordinators at the JRC. … The
> EU-HYDI database is not distributed outside the participating institutions.
> External partners can access only derivatives for joint publications.

So it cannot be downloaded at all, by anyone outside the consortium — a
stricter bar than HYPRES's. The same section does, however, name a legitimate
route:

> scientists from any institution are welcome to contact any author of this
> report for cooperative research. In such research projects handling and
> analysis of raw data has to be done by the EU-HYDI contributing participant
> without giving access to the raw data to the external partner. Results of the
> analysis and derived information can be published together with external
> partners.

One openly licensed European source *was* found and tested: **sDB** (Vereecken,
Van Looy, Weynants & Javaux 2017, doi:10.1594/PANGAEA.879233, **CC-BY-3.0**) —
38 Belgian profiles / 182 horizons, of which 165 yield usable van Genuchten
fits (median 18 retention points per curve, median RMSE 0.014, 136 with Ksat).
Being CC-BY it is redistributable, so `data/sdb_reference.csv` and
`prepare_sdb.py` are included. Measured effect on the 778 European targets,
paired and profile-level: **+0.3 pp class (p = 0.83), -0.4 pp group (p = 0.69)**
— no benefit. The reason is straightforward: GSHP already ingests 145 layers
referenced `Vereecken_et_al_2017`, so ~88 % of sDB is in the reference already.
It is therefore not merged by default.

**Conclusion.** There is no openly licensed, sample-level European retention
database that adds materially to what is already in use, and the open sources that do exist (GSHP, UNSODA) are
already integrated — GSHP in particular already carries the European
contributions of Nemes et al. 2001, Vereecken et al. 2017, Stolbovoy et al.
2016, Richard & Lüscher and Schindler & Müller. Closing the European gap
therefore requires a *collaboration* with an EU-HYDI participant, who would run
this tool's evaluation on their side, rather than any data transfer. Failing
that, the gap is a documented limitation, which is why the tool reports class
probabilities and uncertainty ranges instead of a single answer.

### Obtaining HYPRES (for reference)

Access (ESDAC, JRC): the [European Soil Database v2.0 page]
(https://esdac.jrc.ec.europa.eu/content/european-soil-database-v20-vector-and-attribute-data)
carries an online request form asking for name, e-mail, organisation and type,
country, and a purpose statement (>= 30 characters); download instructions
follow by e-mail. HYPRES is described as a component of the ESDB but is not
offered as a standalone download, so state explicitly that you need the
**sample-level HYPRES tables (RAWRET / RAWK / HYDRAULIC_PROPS)**. The original
custodian is Wageningen Environmental Research (successor to the Winand Staring
Centre; the JRC metadata names J.H.M. Wösten as contact).

Licence, in short: registration is required; ESDB v2.0 terms prohibit passing
data to third parties and prohibit commercial use. **This is incompatible with
redistributing a derived HYPRES table from this public CC BY 4.0 repository**,
so `data/hypres_reference.csv` and `data/hypres*` are git-ignored. Keep the
data local unless you obtain written permission to redistribute a derivative.

Once you have the data:

    python prepare_hypres.py --inspect <file-or-dir>   # see detected columns
    python prepare_hypres.py <file-or-dir>             # write data/hypres_reference.csv

`prepare_hypres.py` accepts .csv/.xlsx/.mdb, uses fitted van Genuchten
parameters when present and otherwise fits the raw theta(h) pairs, and matches
columns by pattern — override anything mis-detected with `--map alpha=ALFA`,
and set `--alpha-units/--head-units/--ksat-units` if they differ from the
defaults (1/cm, cm, cm/day). It prints class-median alpha and n against GSHP as
a unit check: a ~10x offset means the units flag is wrong.

To measure the benefit, add the table as a third variant in
`verify_variants_large.py` (concatenate it exactly as UNSODA is) and re-run the
paired profile-level leave-one-out, ideally reporting European targets
separately — that is where the gain is predicted to land.

### Adding NCSS/KSSL (prepare_kssl.py, verify_kssl.py)

The USDA-NRCS **NCSS Soil Characterization Database** (Kellogg Soil Survey
Laboratory) is public domain and needs no registration. The full SQLite
snapshot (1.95 GB zip → 5.4 GB database) is at
[ncsslabdatamart.sc.egov.usda.gov](https://ncsslabdatamart.sc.egov.usda.gov/database_download.aspx);
`prepare_kssl.py` rebuilds `data/kssl_reference.csv` from it.

**It contains no Ksat.** `lab_analyte` declares "Saturated Conductivity,
Replicate 0/1/2" and "Hydraulic Conductivity", but their `column_name` is empty
and no table in the snapshot holds the values — the same "documented but not
distributed" pattern as HYPRES. Every KSSL row therefore has `ksat_cmh` blank.
Note also that `lab_rosetta_key` must **not** be used as reference van
Genuchten parameters: those are PTF predictions derived *from* texture, so
using them would make the reference circular.

**Size.** Only ~4.3 k of 420 k physical-property layers carry ≥5 retention
tensions plus particle-size analysis — most KSSL pedons store only the
33/1500 kPa pair. After quality control, **2,508 layers** are usable. But the
class mix is complementary to GSHP's: 1.8× the silty clay loam, ~1× the
loam / silt loam / clay loam, and only 4 % as much sand — it fills the classes
where GSHP is starved. (Silt stays starved: 32 + 4.)

**The wet-end truncation artifact.** KSSL's wettest tension is 6 kPa, so
nothing constrains saturation. Fitted as-is, alpha came out 5–16× lower than
GSHP and 19.1 % of layers had θs > 0.75, which is not physical for a mineral
soil. This is *not* a unit error: truncating GSHP's own measured curves at
6 kPa and refitting reproduces the offset in the same direction and magnitude
(clay 0.17×, clay loam 0.17×, silty clay loam 0.20×). The fix is one synthetic
anchor at h = 0 with θ = porosity = 1 − BD_od/2.65, floored just above the
wettest measurement (oven-dry porosity falls *below* it in shrink-swell clays).
Validated on GSHP, where the full-curve answer is known: per-layer agreement
with the true alpha improves from 20 % to 45 % within a factor of two, and
median bias falls from +0.57 to +0.27 dex. After anchoring, θs > 0.75 drops to
0.5 % and alpha aligns with GSHP (loam 0.194 vs 0.197, clay loam 0.120 vs
0.098, silt loam 0.054 vs 0.081). It removes most of the systematic bias but
not the per-layer scatter, so KSSL rows remain noisier than GSHP rows.

**Unit consistency across datasets.** Because α = 1/h_b, the suction unit
scales α directly, so every loader was audited — dimensionally *and*
empirically. Refitting GSHP layers from their raw (h, θ) points reproduces the
published α with a median factor of **1.00** (IQR ±0.02) in every class,
confirming `lab_head_m` is metres and `prepare_gshp.py`'s ÷9.80665 is correct.
KSSL's gravimetric→volumetric conversion was checked independently of α by
comparing θ at 33 and 1500 kPa against GSHP per class (ratios 1.04–1.26, no
BD-sized offset), and its PSDA-derived class matches KSSL's own `texture_lab`
code for **99.6 %** of layers.

| conversion | dataset |
|---|---|
| α[m⁻¹] ÷ 9.80665 | GSHP |
| α[cm⁻¹] ÷ 0.0980665 | ROSETTA |
| fit h in kPa after cm × 0.0980665 | UNSODA, sDB |
| fit h in kPa natively (bar × 100) | KSSL |

**Result** (paired profile-level leave-one-out, 1,663 stratified GSHP targets):

| | GSHP | +KSSL | delta | p |
|---|---|---|---|---|
| exact class | 38.2 % | 39.4 % | +1.2 pp | 0.090 |
| group | 61.5 % | 61.9 % | +0.4 pp | 0.626 |
| Ksat within 2× | 43 % | 44 % | +1 pp | — |
| Ksat RMSE (log10) | 0.90 | 0.89 | −0.01 | — |

The gain is small and **not significant**, though the per-class pattern is
mechanistically coherent — it lands in exactly the classes KSSL enriches
(sandy clay loam +4.0, sandy loam +3.3, clay loam +2.7, silty clay loam +2.7).
The feared Ksat dilution did **not** materialise: even though Ksat-bearing rows
fall from 66 % to 53 % of the reference, ~22 of 30 neighbours still carry a
measured Ksat and every Ksat metric is unchanged.

Exposed as `--reference kssl` (and `all`, which adds UNSODA too). **Defaults
are unchanged**, since neither merge reaches significance on its own.

## Method levers (not more data)

Adding databases gave diminishing returns (UNSODA +0.7 pp p=0.052, KSSL
+1.2 pp p=0.090), so three changes to the *method* were tested instead.

### How optimistic are our numbers? (verify_source_blocked.py)

Every accuracy figure above comes from leave-one-out with the target's
*profile* removed. That rules out site leakage but not **source** leakage:
GSHP is a compilation, and Florida_database alone supplies 57 % of the curated
layers, so a Florida target is still matched against thousands of layers from
the same lab, protocol and region. Real users are not in that position.

Holding out an entire contributing database:

| source | size | n | profile LOO | source-blocked | delta | p |
|---|---|---|---|---|---|---|
| Florida_database | 5735 | 104 | 44.2 % | 27.9 % | −16.3 pp | 0.002 |
| WOSIS | 1011 | 108 | 29.6 % | 20.4 % | −9.3 pp | 0.110 |
| Russia_EGRPR | 1001 | 115 | 20.0 % | 5.2 % | −14.8 pp | 0.002 |
| HYBRAS | 605 | 70 | 32.9 % | 18.6 % | −14.3 pp | 0.013 |
| ETH_Literature | 522 | 99 | 54.5 % | 17.2 % | −37.4 pp | <0.001 |
| AfSPDB | 501 | 111 | 23.4 % | 20.7 % | −2.7 pp | 0.508 |
| UNSODA | 146 | 88 | 40.9 % | 40.9 % | +0.0 pp | 1.000 |
| Belgium_database | 144 | 57 | 35.1 % | 35.1 % | +0.0 pp | 1.000 |
| **POOLED** | 9996 | 925 | **33.8 %** | **22.5 %** | **−11.4 pp** | **<0.001** |

Group accuracy falls 59.4 % → 50.1 % (−9.3 pp, p<0.001). Ksat degrades harder:
within a factor of 2 drops from 25 % to 12 %, Spearman 0.56 → 0.24, and
interval coverage 72 % → 53 %.

**Read the headline numbers accordingly.** Against a soil from an unseen
laboratory, expect roughly **22 % exact class and 50 % group**, not 38/62 %.
Two checks that the design is sound: UNSODA and Belgium show exactly 0.0 pp,
as expected since GSHP ingests near-duplicates of those horizons from other
contributors, so blocking them removes nothing; and AfSPDB, the most
methodologically heterogeneous source, leaks least.

### Tempering the class prior (verify_tau.py, `--tau`)

Votes are weighted by (1/class_count)^tau. At the long-standing tau = 1 the
prior is uniform over the 12 classes, which maximises recall for rare classes
but wrecks their precision: GSHP holds 32 silt layers against 3,930 sand, so
one silt neighbour outvotes 123 sand neighbours.

| tau | overall | macro-R | macro-P | macro-F1 | group | silt precision | silt predicted |
|---|---|---|---|---|---|---|---|
| 0.00 | 38.9 % | 37.0 % | 43.8 % | 34.2 % | 63.1 % | 80.0 % | 5 |
| 0.25 | 39.5 % | 38.4 % | 41.1 % | 36.6 % | 62.9 % | 53.3 % | 15 |
| 0.50 | 39.9 % | 39.7 % | 39.9 % | 38.3 % | 62.9 % | 41.4 % | 29 |
| **0.75** | **39.9 %** | 39.9 % | 39.3 % | **39.0 %** | 62.4 % | 28.3 % | 46 |
| 1.00 (default) | 39.3 % | 40.0 % | 38.7 % | 38.6 % | 61.8 % | 19.5 % | 82 |

There are 32 true silt soils. At tau = 1 the tool answers "silt" 82 times; at
tau = 0.5, 29 times. Macro-recall barely moves (40.0 → 39.7), so the precision
gain is close to free, and Ksat is flat across the sweep. The overall gain is
not significant (+0.7 pp, p=0.33), so **the default stays tau = 1**; use
`--tau 0.75` for best macro-F1 or `--tau 0.5` to calibrate the rare classes.

Note the evaluation is stratified, which makes the target population roughly
uniform — the very prior tau = 1 assumes. This comparison therefore *flatters*
tau = 1; on a natural-frequency evaluation, lowering tau would look better.

### Curve-space matching — tested and rejected (verify_features.py)

The expectation was that measuring distance in (log alpha, log(n−1), thetar,
thetas) distorts curve similarity, because alpha and n trade off along the fit
ridge. Matching on theta at eight fixed potentials (1, 3, 10, 33, 100, 330,
1000, 1500 kPa) should have been closer to a true curve distance.

It is not. Paired against the current default:

| features | tau | overall | macro-R | macro-P | macro-F1 | group | vs default |
|---|---|---|---|---|---|---|---|
| vg (default) | 1.00 | 39.3 % | 40.0 % | 38.7 % | 38.6 % | 61.8 % | — |
| curve | 1.00 | 38.7 % | 38.2 % | 37.9 % | 37.4 % | 62.2 % | −0.6 pp, p=0.60 |
| curve whitened | 1.00 | 37.9 % | 37.7 % | 37.4 % | 36.9 % | 61.6 % | −1.3 pp, p=0.14 |
| vg | 0.75 | 39.9 % | 39.9 % | 39.3 % | 39.0 % | 62.4 % | +0.7 pp, p=0.33 |
| curve whitened | 0.75 | 39.2 % | 38.6 % | 38.0 % | 37.7 % | 62.7 % | −0.1 pp, p=1.00 |

**The premise was backwards.** The eight curve features carry only ~1.6
effective dimensions (PC1 = 77.5 % of variance, mean |correlation| 0.73)
because a retention curve is smooth and every head reports much the same
thing — how wet the soil is. Standardized Euclidean distance over them
collapses to a wetness metric that discards shape, which is why silt, a
shape-defined class, lost 15–19 pp. The vG parameters are the *better*
conditioned space: 3.74 effective dimensions of 4, mean |correlation| 0.12.

Whitening the curve features (`feature_mode="curve_white"`) restores full rank
and fixes the Ksat regression, but makes texture slightly worse still: it
amplifies the near-degenerate directions, which in a smooth curve are mostly
fit noise. Both modes remain available for future experiments; the default is
unchanged.

### Gradient boosting and the accuracy ceiling (verify_gbm.py, verify_ceiling.py)

**Gradient boosting is the one intervention that works.** Compared with kNN on
identical folds and identical features:

| split | model | overall | macro-F1 | group | vs kNN | p |
|---|---|---|---|---|---|---|
| profile-grouped 5-fold | kNN | 36.7 % | 36.0 % | 60.3 % | — | — |
| | kNN + fit-quality wt | 36.7 % | 35.9 % | 60.4 % | +0.0 pp | 1.00 |
| | **GBM** | **39.4 %** | **38.3 %** | 59.8 % | **+2.6 pp** | **0.019** |
| source-blocked | kNN | 17.6 % | 17.0 % | 48.6 % | — | — |
| | **GBM** | **23.0 %** | **19.7 %** | **51.1 %** | **+5.4 pp** | **1.3e-07** |

The gain is largest and most significant exactly where it matters — against a
laboratory the model has never seen, GBM is 31 % better in relative terms. Per
class it wins on clay (+15.3 pp), silt loam (+11.3), silty clay loam (+8.0) and
sandy clay loam (+7.3), and loses on silt (−18.8), silty clay (−12.2) and loam
(−7.3). Its Ksat regressor is slightly better than the neighbour median
(log10 RMSE 0.86 vs 0.91, Spearman 0.68 vs 0.66) but produces **no prediction
interval**, which the kNN does.

**Down-weighting by fit quality does nothing** (+0.0 pp, p=1.00). GSHP
publishes the standard error of its own vG fit, but as a *relative* error the
median is 0.11 and only 0.9 % of layers exceed 1.0 — the reference is mostly
well constrained, so there is little to down-weight. (An earlier reading of the
*absolute* se suggested ~10 % of the reference was unusable; that was wrong,
because absolute se scales with alpha.) Available as `quality_weight=True`,
not recommended.

**The ceiling.** Cover & Hart (1967) bound the Bayes error by the
1-nearest-neighbour error rate, giving a model-free limit for any classifier
built on these four features:

| split | 1NN accuracy | Bayes accuracy bracket | best model today |
|---|---|---|---|
| profile-grouped, exact class | 31.3 % | **[31.3 %, 54.2 %]** | GBM 39.4 % |
| profile-grouped, group | 57.1 % | **[57.1 %, 74.0 %]** | kNN 60.3 % |
| source-blocked, exact class | 19.6 % | **[19.6 %, 40.5 %]** | GBM 23.0 % |
| source-blocked, group | 46.3 % | **[46.3 %, 65.0 %]** | GBM 51.1 % |

Per-class 1NN separability (profile-grouped) ranks sand 83.3 %, clay 52.0 %,
sandy loam 41.3 %, sandy clay loam 34.7 %, silt loam 28.7 %, loam 24.7 %,
loamy sand 20.7 %, silty clay 19.1 %, silty clay loam 18.0 %, clay loam 14.0 %,
sandy clay 11.3 %, silt 9.4 %.

Three caveats. Only the **overall** bracket is rigorous — Cover & Hart bounds
total error, so per-class brackets are heuristic and a few classes fall outside
their own. The bounds are **asymptotic**, so `Bayes >= 1 - err_1NN` holds but
the upper end is an estimate, and the true ceiling may sit a little above it.
And Bayes error depends on the class prior, so both the natural (39 % sand) and
balanced priors are reported.

Note the GBM's own mean top probability (46.5 %) must **not** be read as a
ceiling: its calibration table shows systematic overconfidence (77.6 % mean
confidence in the top bin against 65.9 % actual accuracy).

### The hybrid model (`--model hybrid`, verify_hybrid.py)

Gradient boosting predicts the class better, but a tree ensemble gives no
prediction interval, no particle fractions and no list of similar real soils.
The hybrid keeps both halves: **the GBM supplies the texture class, the kNN
still supplies everything else.**

Measured end to end through `estimate()`, profile-grouped 5-fold:

| model | overall | macro-R | macro-P | macro-F1 | group |
|---|---|---|---|---|---|
| kNN | 36.7 % | 37.4 % | 36.0 % | 36.0 % | 60.3 % |
| **hybrid** | **39.9 %** | **39.0 %** | **39.7 %** | **38.8 %** | 60.1 % |

**+3.2 pp, McNemar p = 0.0047.** That is slightly better than the GBM scored in
isolation (39.4 %), because the shipped hybrid averages `predict_proba` over
the Monte Carlo draws of the fit covariance rather than predicting once from
the point estimate — the same draws the kNN votes over, so fit uncertainty is
still propagated, with a mild ensembling benefit on top.

Per class it gains clay +16.0 pp, silt loam +10.7, silty clay loam +10.7,
sandy clay loam +7.3, and loses silt −18.8, silty clay −9.2, loam −6.0. The
silt loss is less bad than it looks: the kNN's high silt recall came with 19.5 %
precision, and macro-precision rises 36.0 → 39.7 overall.

The verification asserts that **Ksat and its interval are bit-identical**
between the two models, so the GBM demonstrably touches only the class.

**A free confidence signal.** Because the kNN still runs, its own answer is
reported as a second opinion. They agree on 58.9 % of soils, and agreement is
strongly informative:

| | n | hybrid accuracy |
|---|---|---|
| kNN agrees with GBM | 979 | **47.8 %** |
| kNN disagrees | 684 | **28.7 %** |

A 19-point spread, free of charge. Treat a disagreement as a flag that the
answer is unreliable.

Training costs ~3 s on the 9,996-layer reference (early stopping settles near
50 iterations), so the model is fitted on demand rather than shipped as a
pickle, which would tie the repository to one scikit-learn version.

### Where this leaves things

Seven interventions have now been tested on the same paired design. Six were
flat — two data merges (UNSODA +0.7 pp p=0.052, KSSL +1.2 pp p=0.090), three
metric changes (tau +0.7 pp p=0.33, curve −0.6 pp, curve whitened −1.3 pp) and
fit-quality weighting (+0.0 pp p=1.00). One worked: **gradient boosting**,
+2.6 pp in-distribution (p=0.019) and +5.4 pp against an unseen laboratory
(p=1.3e-07).

So the picture is no longer "everything is at the ceiling". It is:

* **The lookup was leaving real accuracy on the table**, and a discriminative
  model recovers a meaningful part of it, especially out-of-lab.
* **But the ceiling is genuinely low.** No classifier on these four van
  Genuchten features can exceed roughly **54 % exact class** in-distribution or
  **40 %** against an unseen laboratory. Group-level ceilings are 74 % and
  65 %. More reference data cannot move those numbers; only more informative
  *inputs* can.
* **The honest headline figure is out-of-lab, not leave-one-out.** Quote
  ~23 % exact class and ~51 % group for a soil from a new source.

The hybrid described above is the result: `--model hybrid` reaches 39.9 %
in-distribution against a ~54 % hard ceiling, keeping every kNN by-product.

Given that ceiling, the highest-value remaining work is about *reporting*
rather than point accuracy: the Ks p5–p95 interval is overconfident (covering
72–84 % of measurements against a nominal 90 %, and only 53 % source-blocked),
and calibrated prediction *sets* would turn irreducible ambiguity into a
defensible output — the kNN/GBM agreement flag is already a crude version of
that, separating 47.8 % from 28.7 % accuracy.

**The default is still `--model knn`**, only because `hybrid` requires
scikit-learn, which the tool otherwise does not need. On accuracy grounds the
hybrid should be the default.

## References

- Soil Survey Staff, Natural Resources Conservation Service, United States
  Department of Agriculture. *National Cooperative Soil Survey Soil
  Characterization Database* (Kellogg Soil Survey Laboratory).
  https://ncsslabdatamart.sc.egov.usda.gov/ (accessed 2026-08-09). Public
  domain (US Government work).
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
- Vereecken, H., Van Looy, K., Weynants, M. and Javaux, M. (2017). Soil
  retention and conductivity curve data base sDB, link to MATLAB files.
  PANGAEA, doi:10.1594/PANGAEA.879233. Licensed **CC-BY-3.0**; redistributed
  here as `data/sDB.zip` and the derived `data/sdb_reference.csv` under that
  licence, with attribution to the authors.
- Nemes, A., Schaap, M.G., Leij, F.J. and Wösten, J.H.M. (2001). Description of
  the unsaturated soil hydraulic database UNSODA version 2.0. *Journal of
  Hydrology* 251:151–162. Data: Ag Data Commons / figshare 24851832.
  (Second reference database, used by `--reference merged`.)
- Wösten, J.H.M., Lilly, A., Nemes, A. and Le Bas, C. (1999). Development and
  use of a database of hydraulic properties of European soils (HYPRES).
  *Geoderma* 90:169–185. Class PTFs obtained from the European Soil Database
  v2.0 (ESDAC, JRC). Only the class and continuous pedotransfer functions were
  ever released, so HYPRES serves here as a European benchmark
  (`verify_hypres.py`), not as reference data.
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
