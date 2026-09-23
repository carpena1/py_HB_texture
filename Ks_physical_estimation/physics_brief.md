# Physics-based estimates from the retention curve: brief

Status on 2026-09-23, commit `1e01b5e`. Code: `ks_physical.py` (the
calculation), `verify_ks_physical.py` (the tests),
`figures/make_external_report.py` (fig6 and `ks_variants.csv` per external
set).

## 1. What it is for

The tool's Ks comes from a lookup, the ~30 nearest reference soils that carry
a measured Ks. The physics gives a second estimate from the curve alone, by
capillary theory, with no neighbours behind it. Its errors are independent
of the reference's, so it serves both as a cross-check and as an answer
where the reference has no Ks (for example `--reference kssl`).

## 2. The calculation

The soil is treated as a bundle of capillaries (Childs and Collis-George
1950; Marshall 1958; Hillel 1980, ch. 8-A, 8-B, 8-I, 9-G).

1. **Suction to pore radius** (capillary equation, contact angle 0):
   r = 2σ / (ρ g h).
2. **Pore-size distribution from the fitted curve.** The range θr–θs is cut
   into m = 100 equal water-content increments. The suction at each
   increment's midpoint is read from the inverse van Genuchten (1980) curve:
   h = [Se^(−1/m′) − 1]^(1/n) / α, with m′ = 1 − 1/n (or the fitted m when m
   is free).
3. **Each pore conducts by Hagen–Poiseuille**, and Marshall's pairing
   statistic sums them, with the radii in decreasing order:

   k = (ε² / 8m²) Σᵢ (2i − 1) rᵢ²,  ε = θs − θr,  Ks = k ρ g / η

   The widest pore carries weight 1 and the narrowest weight 2m − 1.

Constants are for water at 20 °C: σ = 0.0728 N/m, η = 1.002 × 10⁻³ Pa s,
ρ = 998.2 kg/m³, g = 9.81 m/s². The output is in cm/h.

### Correction 1: air-entry cap (the one that matters)

Uncapped, the theory runs high by a factor of six (+0.76 dex). For n < 2 the
van Genuchten curve implies ever wider pores as Se → 1, and the few widest
increments dominate the r² sum. On Zanjanrood (α ≈ 1.3 kPa⁻¹, n ≈ 1.2), for
example, the widest increment has r ≈ 3.6 mm. This is the known near-
saturation defect of the vG form (Vogel et al. 2001; Ippisch et al. 2006).

The fix is to count no pore wider than the one that empties at h_min:
h ← max(h, h_min). The cap was **chosen out of fold**. For each of the 19
scored sources, it was picked from a grid (none, 2, 5, 10, 20, 30, 50, 100
and 200 cm) as the value whose matched physics did best on the other
sources, with every factor refitted without the scored source:

- The folds chose 50 cm in 13 cases, 100 cm in 5 and 200 cm in 1.
- The nested score was 0.61 dex against the kNN's 0.79 (p < 10⁻³⁰).
- **Shipped: h_min = 50 cm** (`AIR_ENTRY_CM`). That is a pore of radius
  30 µm, about 60 µm across, close to the 50 µm lower limit of Greenland's
  (1977) transmission pores.

### Correction 2: matching factor

Classical theory scales the result to a measured Ks, because the curve fixes
the shape of the pore-size distribution but not its absolute level (Childs
and Collis-George 1950; Jackson 1972). Here the factor is:

- **Definition:** the median, over laboratories, of each laboratory's median
  ratio of physics to measured Ks. Each source gets one vote, so Florida,
  which holds half the layers, does not set it.
- **With the cap:** 1.085 on the merged reference (23 sources with a
  measured Ks) and 1.087 on the public one (13). **Shipped: 1.09**
  (`MATCHING_FACTOR`). The laboratories' own ratios run from 0.15 to 200.
- **Without the cap:** 3.9, with the laboratories running from 0.8 to 8,000.
  Weighted by layer instead of by source it would be 2.4.

With the cap, the theory's level is nearly right with no scaling at all.

### Alternative tested: Peters et al. (2023)

Peters et al. write the bundle in closed form so that no matching factor is
needed: Ks = β τs (θs − θr)² α², with β = 3.04 × 10⁻⁴ m³/s and τs = 0.062
for the constrained vG curve. It is implemented as `peters_ks`, and it
scored worse than capped Marshall (Section 3).

## 3. Validation

**Source-blocked.** These are 1,870 targets from the 19 sources with at
least 40 measured-Ks layers, stratified by class. Each source was hidden from
the kNN and the classifier, curves were regenerated from the stored vG
parameters and refitted, and the factor was refitted without the source.

| Ks estimate | Median error | Within 2× / 10× | Bias | ρ |
|---|---|---|---|---|
| kNN (the tool's Ks) | 0.79 dex | 23 % / 59 % | +0.09 | 0.29 |
| Physics, 50 cm cap, raw | 0.62 | 28 % / 69 % | +0.07 | 0.51 |
| **Physics, 50 cm cap, matched (shipped)** | **0.61** | 28 % / 68 % | +0.04 | 0.51 |
| Physics, 10 cm cap, matched | 0.68 | 25 % / 65 % | +0.07 | 0.47 |
| Physics, no cap, raw | 1.00 | 16 % / 50 % | +0.76 | 0.44 |
| Physics, no cap, matched | 0.78 | 21 % / 59 % | +0.16 | 0.43 |
| Peters et al. 2023 | 0.91 | 18 % / 54 % | +0.40 | 0.39 |

The ρ column is the Spearman rank correlation with the measured Ks.

- **By source:** the shipped physics beats the kNN in 16 of the 19 sources.
  The exceptions are the Yellow River (+0.04 dex), EU-HYDI Romano (+0.02) and
  EU-HYDI Lilly (+0.36; its θ(0) is 1.07 × porosity).
- **By texture group** (kNN → physics, dex): sandy 0.45 → 0.34, loamy
  0.75 → 0.63, silty 1.09 → 0.80, clayey 0.98 → 0.68.
- **Agreement as a confidence signal.** The matched physics falls within 10×
  of the kNN for 79 % of soils. There the kNN's median error is 0.65 dex,
  against 1.68 where they disagree (67 % vs 29 % within 10×).
- **Blending** adds nothing. With the weight fitted on the other sources it
  goes to 0.95 on the physics and scores 0.61 dex, and a fixed half-and-half
  is worse (0.67).

**External sets** (typical error factor; the set is removed from the
reference):

| Set | n | kNN | Physics, shipped | Physics, no cap, matched |
|---|---|---|---|---|
| Zanjanrood (Iran) | 169 | ×2.35 | **×1.95** | ×16.8 |
| Arizona | 21 | ×2.0 | ×2.0 | ×2.5 |
| Laikipia (field permeameter Ks) | 84 | ×11.5 | ×10.2 | ×11.6 |
| New Jersey, SSIR 26 (intact-core Ks) | 236 | ×1.84 | ×2.22 (ρ 0.80) | ×4.99 |

The New Jersey Ks played no part in choosing the cap or the factor: uncapped and raw, the physics would be ×16.9 there. Boorowa has no measured Ks. Laikipia's Ks is a field measurement that no
lab-based estimate reproduces, so it cannot separate the methods.

**A range for the physics was tested and rejected.** The spread across the
variants holds the measured Ks only 10–39 % of the time. The honest
empirical 5–95 % range is 3.7 dex wide (×/÷ 70) and the same for every soil.

## 4. How the tool reports it (user decisions, 2026-09-23)

```
Predicted Ks: 2.01 cm/h  [5-95 %: 8.33e-05 - 34.7]   (from the kNN; ~12 of 30 undisturbed neighbors with measured Ksat)
  physical second opinion agrees: 0.42 cm/h, 4.8x below  (0.458 raw; Marshall 1958 capillary bundle, pores capped at 50 cm suction, matched = raw / 1.09)
```

- **The kNN stays the Ks, with its 5–95 % band.** It is the only estimate
  with a band, and it improves when the user's own laboratory is in the
  reference (0.53 dex for a new site from a known source, 0.40 for a new
  depth). The physics cannot use that.
- **The physics is printed beside it, raw and matched, with no band.**
  "Agrees" means within a factor of 10 of the kNN, judged on the matched
  value.
- **Not adopted:** the physics median with the kNN band, blending, a
  variant range, and the Peters form.

## 5. Texture from the pore sizes (tested, not adopted)

The bundle can also be run toward texture: invert Arya and Paris (1981),
pore radius → particle radius with their scaling exponent α = 1.38, then
cumulative particle-size distribution → USDA fractions and class. There is
no fitting. It was tested paired against the tool on 1,733 balanced,
source-blocked targets:

| | Exact | Group | Macro-F1 | Fraction error (points) |
|---|---|---|---|---|
| Tool (classifier) | 31.2 % | 57.9 % | 29.9 | 13.3 |
| Neighbour vote | 28.2 % | 55.5 % | 27.9 | – |
| Inverse Arya–Paris | 21.2 % | 45.0 % | 17.1 | 18.2 |

It is −10.0 pp against the tool (p = 2 × 10⁻¹⁶). It nearly matches the tool
on sandy soils (58 % vs 61 %) and collapses on clays (10 % vs 27 %; clay
content error 15.9 vs 8.5 points), where aggregation and film flow break the
capillary assumption. An earlier figure of 35.1 % came from the reference's
sand-heavy natural class mix and was withdrawn.

## 6. Caveats

- **The cap and the factor were chosen on the reference's measured Ks.**
  The cap was chosen out of fold, but only the "raw" value is free of any
  scaling to measured Ks.
- **Zanjanrood's Ks is from repacked samples.** It correlates with nothing
  in the table (|ρ| < 0.1), so it tests the Ks level only. Its intact-core
  wet end drains 54 % of its water range by 100 cm (26 % for other
  undisturbed clay loams), and the curves are spliced from two samples at
  100/330 cm.
- **EU-HYDI Kätterer** is offset from every other source for both
  estimates (kNN 2.54 dex, physics 2.12). That looks like a method or unit
  difference in the source.
- **The theory is weakest in fine soils, and so is the kNN.** Capillary
  flow ignores aggregation, film flow and tortuosity differences, but with
  the cap the physics still beats the kNN in every texture group.

## 7. Reproduce

```bash
.venv/bin/python verify_ks_physical.py                 # ~10 min: all tables above
.venv/bin/python figures/make_external_report.py babaeian_zanjanrood   # fig6, ks_variants.csv
```

Wrap long runs in `caffeinate -i`, with OMP/OPENBLAS/VECLIB/MKL_NUM_THREADS=1.

## References

- Arya, L.M. and Paris, J.F. (1981). A physicoempirical model to predict the
  soil moisture characteristic from particle-size distribution and bulk
  density data. *Soil Science Society of America Journal* 45:1023–1030.
- Childs, E.C. and Collis-George, N. (1950). The permeability of porous
  materials. *Proceedings of the Royal Society of London A* 201:392–405.
- Greenland, D.J. (1977). Soil damage by intensive arable cultivation:
  temporary or permanent? *Philosophical Transactions of the Royal Society
  of London B* 281:193–208.
- Hillel, D. (1980). *Fundamentals of Soil Physics*. Academic Press, New
  York. Chapters 8 and 9.
- Ippisch, O., Vogel, H.-J. and Bastian, P. (2006). Validity limits for the
  van Genuchten–Mualem model and implications for parameter estimation and
  numerical simulation. *Advances in Water Resources* 29:1780–1789.
- Jackson, R.D. (1972). On the calculation of hydraulic conductivity. *Soil
  Science Society of America Proceedings* 36:380–382.
- Marshall, T.J. (1958). A relation between permeability and size
  distribution of pores. *Journal of Soil Science* 9:1–8.
- Peters, A., Hohenbrink, T.L., Iden, S.C., van Genuchten, M.Th. and
  Durner, W. (2023). Prediction of the absolute hydraulic conductivity
  function from soil water retention data. *Hydrology and Earth System
  Sciences* 27:1565–1582.
- van Genuchten, M.Th. (1980). A closed-form equation for predicting the
  hydraulic conductivity of unsaturated soils. *Soil Science Society of
  America Journal* 44:892–898.
- Vogel, T., van Genuchten, M.Th. and Císlerová, M. (2001). Effect of the
  shape of the soil hydraulic functions near saturation on variably-saturated
  flow predictions. *Advances in Water Resources* 24:133–144.
