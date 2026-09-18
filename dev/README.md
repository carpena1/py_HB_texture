# Development history

Experiments that shaped the current tool but are no longer part of its
verification. Each one tested a change and either led to the current defaults
or was rejected; the main README summarises the outcomes. The results quoted in
their docstrings were computed on the reference and defaults of their day, so
re-running them now will not reproduce those numbers exactly.

Run them from the repository root with the root on the import path:

    PYTHONPATH=. .venv/bin/python dev/verify_tau.py

| script | question it answered |
|---|---|
| `verify_variants.py`, `verify_variants_large.py` | Do UNSODA or depth help? (paired leave-one-out) |
| `verify_by_class.py` | Per-class and per-group accuracy, with and without depth |
| `verify_region.py` | What is region-matched reference data worth? |
| `verify_skill.py`, `verify_skill_hybrid.py`, `verify_skill_kssl.py` | Is each class acceptable in absolute terms? |
| `verify_kssl.py` | Does merging NCSS/KSSL help? |
| `verify_tau.py` | Tuning the class-prior exponent |
| `verify_features.py` | vG parameters vs curve-space matching |
| `verify_gbm.py` | Gradient boosting vs kNN on matched folds |
| `verify_om.py` | Does organic carbon help? |
| `verify_free_m.py` | Does freeing m from n help? |
| `verify_hohenbrink.py` | What the Hohenbrink et al. (2023) cores add; Ks on structured soils |
| `verify_euhydi.py` | What EU-HYDI adds (needs the restricted table) |
