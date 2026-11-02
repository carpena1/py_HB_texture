# Validation scripts

The tests behind the accuracy figures in the main README, the admission test
for new data, and the tests that set the tool's design. Run every script from
the repository root, for example:

    .venv/bin/python validation/verify_holdout.py

Each script states its arguments, its result and its date in its docstring.
The scripts use all cores (`verify_common.pool`, one process per core and one
thread per process); set `SWCC_JOBS` to use fewer. Slow intermediate results
are cached in `figures/cache/`.

## Accuracy at each validation level

| script | question |
|---|---|
| `verify_holdout.py` | How much does it matter what is held out when a test soil is hidden? (the headline accuracy) |
| `verify_source_blocked.py` | How much accuracy survives a soil from a laboratory the reference does not hold? |
| `verify_ceiling.py` | A model-free estimate of the accuracy ceiling |
| `verify_groups.py` | Accuracy for the aggregated texture groups |
| `verify_hybrid.py` | The hybrid model (classifier plus neighbours) as the tool runs it |
| `verify_fractions.py` | Neighbour-mean particle fractions against sand and clay regressors |
| `verify_gshp.py` | Leave-one-out check against real GSHP soils |
| `verify_carsel_parrish.py`, `verify_rosetta.py` | Checks against the Carsel & Parrish and ROSETTA class means |
| `verify_predictive.py` | Is the tool predictive for a new sensor network, and how accurate is it? |
| `verify_network.py` | What a new sensor network can expect, and what improves it |

## New data

| script | question |
|---|---|
| `verify_external.py` | How does the tool read an outside dataset before it joins the reference? |
| `verify_join_effect.py` | Does adding a candidate table make the predictions for the other soils worse? (the admission test) |

## Saturated conductivity (Ks)

| script | question |
|---|---|
| `verify_ks_neighbors.py` | How many neighbours should Ks come from, and of which layers? |
| `verify_ks_physical.py` | Is the capillary-bundle Ks worth reporting beside the neighbours, and how should it be corrected? |
| `verify_literature_methods.py` | Do the methods of the literature review do better, for Ks and for clay? |

## Design tests

| script | question |
|---|---|
| `verify_covariates.py` | Do sample depth and bulk density improve the prediction? |
| `verify_sample_type.py` | Does knowing how a sample was prepared help? |
| `verify_kssl_label.py` | Does KSSL's sample-type label change anything? |
| `verify_alt_predictors.py`, `verify_ivg.py` | Do alternative predictors (IvG parameters, fixed-head water contents) do better? |
| `verify_blend.py` | Can the classifier and the neighbour vote be combined by class? |
| `verify_partial_heads.py` | Can soils with only two or three retention points still teach the tool? |
| `verify_volcanic.py`, `verify_andic.py`, `verify_andic_knn.py` | Volcanic and andic soils |
| `verify_psd_standard.py`, `verify_sieve_1mm.py` | Particle sizes measured on a different sieving basis |

## Shared helpers

- `verify_common.py` — shared constants and statistics, thread pinning and the process pool (`pool`, `merge`, `map_folds`)
- `verify_groups.py` — the texture groups
- `ksat_metrics.py` — scoring for Ks (log10 errors, factor-of-N hits)

## `archive/`

Superseded experiments, kept for the record; see `archive/README.md`.
