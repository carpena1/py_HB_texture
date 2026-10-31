# Physics-based estimates from the retention curve: brief

This brief covers the estimates that use physics (or a published equation)
on the retention curve alone, with no reference soils behind them. Part A is
about saturated hydraulic conductivity (Ks). The tool reports one such
estimate, the capillary bundle, as a second opinion. Part B is about
texture. No physical texture method was good enough, so the tool reports
none; Part B records what was tried.

**Status.** First written on 2026-09-23 (commit `1e01b5e`). Texture methods
updated 2026-10-01. kNN Ks figures updated 2026-10-03 for the new Ks lookup
(150 Ks-bearing neighbours, `verify_ks_neighbors.py`). Methods from the
literature review (`Literature_metanalysis/review.md`) added 2026-10-30
(`verify_literature_methods.py`). CPA on three independent sets added
2026-10-31.

**Code.** `ks_physical.py` (the calculation), `verify_ks_physical.py` (the
Ks tests), `verify_literature_methods.py` (the literature forms, Ks and
clay), `figures/make_external_report.py` (fig6 and `ks_variants.csv` for
each external set).

## Contents

- [Part A. Ks prediction](#part-a-ks-prediction)
  - [A1. What it is for](#a1-what-it-is-for)
  - [A2. The calculation](#a2-the-calculation)
    - [Correction 1: air-entry cap](#correction-1-air-entry-cap)
    - [Correction 2: matching factor](#correction-2-matching-factor)
    - [Alternatives tested: Peters et al. (2023) and Brutsaert (1967)](#alternatives-tested-peters-et-al-2023-and-brutsaert-1967)
  - [A3. Validation](#a3-validation)
  - [A4. Other Ks forms from the literature review](#a4-other-ks-forms-from-the-literature-review)
    - [CPA on three independent sets](#cpa-on-three-independent-sets)
  - [A5. How the tool reports Ks](#a5-how-the-tool-reports-ks)
  - [A6. Caveats](#a6-caveats)
- [Part B. Texture prediction](#part-b-texture-prediction)
  - [B1. What it is for, and the result](#b1-what-it-is-for-and-the-result)
  - [B2. Four reference-free methods](#b2-four-reference-free-methods)
  - [B3. Clay from the dry end](#b3-clay-from-the-dry-end)
  - [B4. Not tested, and why](#b4-not-tested-and-why)
- [Reproduce](#reproduce)
- [References](#references)

---

## Part A. Ks prediction

### A1. What it is for

The tool's Ks comes from a lookup. It takes the 150 nearest reference soils
that carry a measured Ks (`swcc_texture.KS_K`). Until 2026-10 it took the 30
nearest soils of any kind, and only those with a Ks voted.

The physics gives a second estimate from the curve alone. It uses capillary
theory and no neighbours. Its errors do not depend on the reference, so it
has two uses:

- a cross-check on the lookup;
- an answer where the reference has no Ks (for example `--reference kssl`).

### A2. The calculation

The model treats the soil as a bundle of capillaries (Childs and
Collis-George 1950; Marshall 1958; Hillel 1980, ch. 8-A, 8-B, 8-I, 9-G).

1. **Convert suction to pore radius.** Use the capillary equation with
   contact angle 0: r = 2σ / (ρ g h).
2. **Get the pore-size distribution from the fitted curve.** Cut the range
   θr–θs into m = 100 equal water-content steps. Read the suction at the
   midpoint of each step from the inverse van Genuchten (1980) curve:
   h = [Se^(−1/m′) − 1]^(1/n) / α, with m′ = 1 − 1/n (or the fitted m when
   m is free).
3. **Let each pore conduct by Hagen–Poiseuille.** Marshall's pairing
   statistic sums the pores, with the radii in decreasing order:

   k = (ε² / 8m²) Σᵢ (2i − 1) rᵢ²,  ε = θs − θr,  Ks = k ρ g / η

   The widest pore has weight 1. The narrowest has weight 2m − 1.

The constants are for water at 20 °C: σ = 0.0728 N/m,
η = 1.002 × 10⁻³ Pa s, ρ = 998.2 kg/m³, g = 9.81 m/s². The output is in
cm/h.

#### Correction 1: air-entry cap

This is the correction that matters. Without it, the theory is too high by
a factor of six (+0.76 dex). For n < 2, the van Genuchten curve implies
wider and wider pores as Se → 1, and the few widest steps control the r²
sum. On Zanjanrood (α ≈ 1.3 kPa⁻¹, n ≈ 1.2), for example, the widest step
has r ≈ 3.6 mm. This is the known defect of the vG form near saturation
(Vogel et al. 2001; Ippisch et al. 2006).

The fix counts no pore wider than the pore that empties at h_min:
h ← max(h, h_min). The cap was **chosen out of fold**. For each of the 20
scored sources, the value came from a grid (none, 2, 5, 10, 20, 30, 50, 100
and 200 cm). The chosen value was the one whose matched physics did best on
the other sources. Every factor was refitted without the scored source.

- The folds chose 50 cm in 17 cases and 100 cm in 3. (In the first run, on
  19 sources, they chose 50 cm 13 times, 100 cm 5 times and 200 cm once.)
- The nested score was 0.62 dex, against 0.68 for the kNN (p = 10⁻⁷).
  Against the earlier 30-neighbour kNN it was 0.61 against 0.79
  (p < 10⁻³⁰).
- **Shipped: h_min = 50 cm** (`AIR_ENTRY_CM`). That is a pore of radius
  30 µm, about 60 µm across. This is close to the 50 µm lower limit of
  Greenland's (1977) transmission pores.

#### Correction 2: matching factor

The curve gives the shape of the pore-size distribution, but not its
absolute level. So classical theory scales the result to a measured Ks
(Childs and Collis-George 1950; Jackson 1972). Here the factor is:

- **Definition:** the median, over laboratories, of each laboratory's
  median ratio of physics to measured Ks. Each source gets one vote. So
  Florida, which holds half the layers, does not set the factor.
- **With the cap:** 1.085 on the merged reference (23 sources with a
  measured Ks) and 1.087 on the public one (13). **Shipped: 1.09**
  (`MATCHING_FACTOR`). The ratios of single laboratories run from 0.15 to
  200.
- **Without the cap:** 3.9, with laboratories from 0.8 to 8,000. Weighted
  by layer instead of by source, it would be 2.4.

With the cap, the level of the theory is almost right with no scaling.

#### Alternatives tested: Peters et al. (2023) and Brutsaert (1967)

**Peters et al. (2023)** write the capillary bundle in closed form, so that
no matching factor is necessary: Ks = β τs (θs − θr)² α², with
β = 3.04 × 10⁻⁴ m³/s and τs = 0.062 for the constrained vG curve. It is
implemented as `peters_ks`. It scored worse than capped Marshall (A3).

**Brutsaert (1967), with the Brooks–Corey curve.** Brakensiek et al. (1981,
eq. 5) compute Ks from fitted Brooks and Corey (1964) parameters with
Brutsaert's (1967) equation:

   Ks = 270 φe² λ² / (ψb² (λ + 1)(λ + 2))   (cm/s; ψb in cm; φe = φ − θr)

This is not a separate theory. It is the same Childs–Collis-George/Marshall
capillary bundle, in closed form for the Brooks–Corey curve:

- **The constant 270** is σ²/(2μρg) for water at 20 °C in cgs units
  (72.8² / (2 × 0.01 × 1 × 981) = 270 cm³/s). Marshall's sum has the same
  prefactor when each radius is written as r = 2σ/(ρgh).
- **The rest:** with many steps, Marshall's pairing sum becomes
  ∫₀¹ 2(1 − S) / h(S)² dS. With h = ψb S^(−1/λ), this integrates exactly to
  λ² / (ψb² (λ + 1)(λ + 2)).

The one real difference from our van Genuchten form is the curve. The
Brooks–Corey curve is flat above its bubbling pressure ψb, so no pore is
wider than r(ψb). The curve thus has an air-entry cap built in. Van
Genuchten with n < 2 needs one imposed (Correction 1).

Brakensiek et al. then halve Ks to get the Green–Ampt conductivity. They
follow Bouwer's (1966) re-wet conductivity. That halving is for
infiltration models, not a correction to Ks.

The test fitted Brooks–Corey to the same curves (`brutsaert_ks` in
`verify_ks_physical.py`; the tool does not use it). The fitted ψb have a
median of about 13 cm. So Eq. 5 behaves like the van Genuchten bundle with
a cap of a few centimetres. It is better than the form with no cap, but not
as good as the 50 cm cap, both source-blocked and on the New Jersey intact
cores (A3). It is not adopted, but it gives the cap a classical anchor.

### A3. Validation

**Source-blocked.** There are 1,906 targets from the 20 sources with at
least 40 measured-Ks layers, stratified by class. For each target:

1. Hide the target's source from the kNN and the classifier.
2. Regenerate the curve from the stored vG parameters and refit it.
3. Refit the factor without the source.

| Ks estimate | Median error | Within 2× / 10× | Bias | ρ |
|---|---|---|---|---|
| kNN (the tool's Ks) | 0.68 dex | 25 % / 64 % | +0.15 | 0.37 |
| Physics, 50 cm cap, raw | 0.62 | 28 % / 69 % | +0.08 | 0.51 |
| **Physics, 50 cm cap, matched (shipped)** | **0.61** | 28 % / 69 % | +0.04 | 0.51 |
| Physics, 10 cm cap, matched | 0.68 | 25 % / 66 % | +0.04 | 0.48 |
| Physics, no cap, raw | 0.99 | 16 % / 50 % | +0.75 | 0.44 |
| Physics, no cap, matched | 0.77 | 22 % / 60 % | +0.15 | 0.44 |
| Peters et al. 2023 | 0.90 | 18 % / 54 % | +0.39 | 0.39 |
| Brutsaert 1967 (Brooks–Corey), raw | 0.75 | 23 % / 61 % | +0.26 | 0.46 |
| Brutsaert 1967, matched (factor 1.57) | 0.72 | 24 % / 63 % | +0.07 | 0.46 |
| Brutsaert 1967, halved (Bouwer re-wet K) | 0.71 | 25 % / 63 % | −0.04 | 0.46 |

ρ is the Spearman rank correlation with the measured Ks.

- **By source:** the shipped physics beats the kNN in 16 of the 20
  sources. The exceptions are the Yellow River (+0.09 dex), Zanjanrood
  (+0.05), EU-HYDI Romano (+0.10) and EU-HYDI Lilly (+0.54; its θ(0) is
  1.07 × porosity).
- **By texture group** (kNN → physics, dex): sandy 0.42 → 0.34, loamy
  0.64 → 0.62, silty 0.89 → 0.79, clayey 0.91 → 0.68.
- **Agreement is a confidence signal.** The matched physics is within 10×
  of the kNN for 87 % of soils. There, the kNN's median error is 0.60 dex.
  Where they disagree, it is 1.49 dex (69 % against 32 % within 10×).
- **Blending adds nothing that matters.** With the weight fitted on the
  other sources, the weight goes to 0.80 on the physics and the blend
  scores 0.60 dex. That is level with the physics alone (0.61). A fixed
  half-and-half blend scores the same.

**External sets.** The table gives the typical error factor. Each set is
removed from the reference.

| Set | n | kNN | Physics, shipped | Physics, no cap, matched |
|---|---|---|---|---|
| Zanjanrood (Iran) | 169 | **×1.81** | ×1.96 | ×16.7 |
| Arizona | 21 | ×2.05 | ×1.98 | ×2.6 |
| Laikipia (field permeameter Ks)* | 84 | ×11.5 | ×10.2 | ×11.6 |
| New Jersey, SSIR 26 (intact-core Ks) | 236 | ×1.83 | ×2.22 (ρ 0.80) | ×4.83 |

\* Laikipia was measured with the earlier 30-neighbour kNN and not re-run,
because its table was then under revision.

- The New Jersey Ks had no part in the choice of the cap or the factor.
  With no cap and no factor, the physics would be ×16.9 there.
- Brutsaert's Eq. 5, with Brooks–Corey fitted to the measured points of
  each horizon (median ψb 12.7 cm), gives ×5.3 there (ρ 0.63), and ×2.8
  halved.
- Boorowa has no measured Ks.
- Laikipia's Ks is a field measurement. No laboratory-based estimate
  reproduces it, so it cannot separate the methods.

**A range for the physics was tested and rejected.** The spread across the
variants holds the measured Ks only 10–39 % of the time. The honest
empirical 5–95 % range is 3.7 dex wide (×/÷ 70), and it is the same for
every soil.

### A4. Other Ks forms from the literature review

The review (`Literature_metanalysis/review.md`) ranks critical-path
analysis (CPA) first among the Ks models that use the curve: RMSLE 0.87 on
313 Kansas Mesonet cores (Ghanbarian and Patrignani 2026). CPA and four
other forms were scored with `verify_literature_methods.py`, exactly as
that paper writes them:

- **CPA:** Ks = (ρg/μ) dc² / (C F), with C = 53.5 (Skaggs 2011).
  - dc is the pore-throat diameter (d = 4σ/ρgh) at the inflection point of
    the van Genuchten curve. The inflection is read on log h,
    h = m^(−1/n)/α. A variant reads it on linear h, h = m^(1/n)/α.
  - F is the formation factor from porosity (Ghanbarian et al. 2014):
    F = (1 − φc)(φx − φc)/(φ − φc)² below the crossover φx = 0.75, and
    (1 − φc)/(φ − φc) above it. The critical porosity is φc = θr. With
    φc = 0 and φx = 1, this is Archie's F = φ⁻².
- **JKS** (Johnson et al. 1986): Ks = (ρg/μ) Λ²/(8F), with Λ = 0.21 dc.
- **Guarracino (2007):** Ks = 4.65 × 10⁴ φ α² cm/d, with α in cm⁻¹.
- **Mishra and Parker (1990):** Ks = 108 (θs − θr)^2.5 α² cm/s.
- **Effective-porosity Kozeny–Carman** (Ahuja et al. 1989):
  Ks = 764.5 φe^3.29 cm/h, with φe = φ − θ(33 kPa).

How they were scored:

- φ is θs. Where bulk density is known, 1 − ρb/2.65 is a variant.
- Each form is also matched with its own factor, defined as in A2 and
  refitted without the scored source.
- The forms use the stored van Genuchten parameters.
- The targets are the 1,843 source-blocked targets of
  `verify_ks_neighbors.py` (19 sources). So the kNN row is that script's
  150-neighbour value.
- "Sources" counts the sources where the form beats the shipped physics.

| Ks estimate | Median error | RMSLE | Within 2× / 10× | Bias | ρ | Sources |
|---|---|---|---|---|---|---|
| **Physics, 50 cm cap (shipped)** | **0.62** | **1.22** | 27 % / 68 % | +0.04 | 0.51 | – |
| kNN (the tool's Ks) | 0.71 | 1.31 | 24 % / 64 % | +0.15 | 0.35 | 4/19 |
| Kozeny–Carman (Ahuja), φ from ρb, raw | 0.66 | 1.99 | 25 % / 65 % | −0.34 | 0.51 | 6/19 |
| Kozeny–Carman (Ahuja), φ = θs, raw | 0.69 | 1.82 | 25 % / 64 % | −0.43 | 0.50 | 5/19 |
| CPA, log-h inflection, matched (factor 3.0) | 0.78 | 1.52 | 21 % / 60 % | +0.08 | 0.44 | 2/19 |
| CPA, log-h inflection, raw | 0.94 | 1.61 | 18 % / 52 % | +0.55 | 0.44 | 2/19 |
| JKS, raw | 0.78 | 1.52 | 23 % / 60 % | +0.02 | 0.44 | 3/19 |
| Guarracino 2007, matched | 0.91 | 1.67 | 17 % / 54 % | +0.17 | 0.36 | 3/19 |
| Mishra and Parker 1990, matched | 0.93 | 1.71 | 17 % / 53 % | +0.15 | 0.37 | 2/19 |
| CPA, linear-h inflection, raw | 2.28 | 3.18 | 4 % / 19 % | +2.48 | 0.29 | 0/19 |
| Half shipped, half Kozeny–Carman (φ from ρb) | 0.62 | 1.45 | 27 % / 68 % | −0.15 | 0.52 | 9/19 |

None is better than the shipped physics. Every row is significantly worse
on median error (p ≤ 10⁻⁴, paired), except the last, which is level
(p = 0.45).

- **On our data, CPA ranks second, not first.** Its RMSLE is 1.52 here,
  against 0.87 on the Kansas cores; the shipped bundle's is 1.22. CPA reads
  one pore size, at the inflection point. The bundle sums over every pore
  up to the 50 cm cap. On a global reference with many laboratories, the
  whole pore-size distribution carries more information than one point. CPA
  loses in every texture group (all 11,100 layers with a Ks, matched: sandy
  0.21 against 0.18 dex, loamy 0.69 against 0.55, silty 0.90 against 0.69,
  clayey 0.95 against 0.72).
- **Read the inflection on log h.** On linear h the critical throat is
  about 9 times wider (median), and Ks comes out about 300 times too high
  (+2.48 dex).
- **JKS is CPA with another constant.** Λ²/8 = 0.0055 d², against
  d²/53.5 = 0.019 d². Its raw form lands where matched CPA does.
- **Kozeny–Carman is close on the median and worse in the tails.** It needs
  θ(33 kPa) and a porosity, so it is the only candidate that uses more than
  the shape of the curve. Its RMSLE is the worst of the plausible forms. It
  under-predicts by 0.3–0.4 dex. A 50/50 blend with the shipped physics
  adds nothing.
- **Guarracino and Mishra–Parker** put all of Ks into α², as Peters et al.
  (2023) do, and they score close to Peters (0.90 dex in A3). Their levels
  are far off: matched, Guarracino is multiplied by 6.7 and Mishra–Parker
  divided by 6.5. Matched, the order CPA, JKS, Guarracino, Mishra–Parker is
  the same as in the Kansas paper.

#### CPA on three independent sets

Is the gap between the Kansas paper and our reference in the method or in
the data? To find out, the same equations (log-h inflection, φ = θs) were
scored on the fitted curves and measured Ks of three independent sets
(2026-10-31). "Matched" divides by the global factor of the default
reference (×2.7).

| Set | n | Shipped bundle | CPA, raw | CPA, matched | kNN, new source |
|---|---|---|---|---|---|
| Kansas Mesonet (clean Ks)* | 251 | 0.47 | **0.42** (p = 0.32) | 0.53 | 0.65† |
| Kansas Mesonet (every positive Ks) | 313 | 0.50 | **0.43** (p = 0.015) | – | – |
| Fujii et al., Japan and SE Asia | 51 | 0.77 | 0.70 (p = 0.37) | **0.58** | 1.24 |
| Ksat-SSIR-DB, Brazil | 95 | **0.57** | 2.08 | 1.69 | 0.60 |

Median |log10 error| in dex. p is the paired Wilcoxon test against the
shipped bundle. \* Positive Ks with no laboratory note (gaps, cracks,
roots). † On the 126 cores that also have a texture.

- **On the Kansas cores, CPA works as the paper says.** Raw, it has almost
  no bias, and it is level with or slightly better than the bundle. So the
  equations are right, and our reference is different.
- **On the Japanese and Southeast Asian forest soils, CPA is again
  close** (bias +0.15 dex). The bundle under-predicts by 0.74 dex. Both
  physics forms beat the kNN, which is about 16 times too low here (median
  measured Ks 23 cm/h).
- **On the Brazilian soils, CPA fails.** It is about 120 times too high
  (+2.08 dex), within 10× for only 14 % of rows. These soils are mostly
  clays: Latossolos (Oxisols, 41 of 95), Cambissolos (28) and Argissolos
  (21). A likely cause, not tested: in these strongly aggregated soils the
  inflection of the curve sits at the pores between aggregates. Those
  pores are wider than the throat that limits the flow. The bundle stays
  good there (0.57 dex, 75 % within 10×, ρ 0.59).
- **So CPA depends on the site.** It does well on soils like those it was
  built on. It does badly where the inflection of the curve is not the
  critical throat. Over the 19 sources of the global test, the bundle is
  the safer choice.

Not adopted. The shipped capillary bundle stays the second opinion.

### A5. How the tool reports Ks

These are user decisions of 2026-09-23.

```
Predicted Ks: 2.16 cm/h  [5-95 %: 0.00433 - 41.7]   (from the kNN: the 150 nearest undisturbed neighbors with a measured Ksat)
  physical second opinion agrees: 0.42 cm/h, 5.2x below  (0.458 raw; Marshall 1958 capillary bundle, pores capped at 50 cm suction, matched = raw / 1.09)
```

- **The kNN stays the Ks, with its 5–95 % band.** It is the only estimate
  with a band. It also improves when the user's own laboratory is in the
  reference (0.44 dex for a new site from a known source, 0.38 for a new
  depth). The physics cannot use that.
- **The physics is printed beside it, raw and matched, with no band.**
  "Agrees" means within a factor of 10 of the kNN, judged on the matched
  value.
- **Not adopted:** the physics median with the kNN band, blending, a range
  from the variants, the Peters form, Brutsaert's Brooks–Corey form, and
  the forms from the literature review (CPA, JKS, Guarracino,
  Mishra–Parker, effective-porosity Kozeny–Carman; A4).

### A6. Caveats

- **The cap and the factor were chosen on the measured Ks of the
  reference.** The cap was chosen out of fold. But only the "raw" value has
  no scaling to measured Ks.
- **Zanjanrood's Ks is from repacked samples.** It correlates with nothing
  in the table (|ρ| < 0.1), so it tests the Ks level only. Its intact-core
  wet end drains 54 % of its water range by 100 cm (26 % for other
  undisturbed clay loams). Its curves join two samples at 100/330 cm.
- **EU-HYDI Kätterer** is offset from every other source for both
  estimates (kNN 2.52 dex, physics 2.11). This looks like a method or unit
  difference in the source.
- **The theory is weakest in fine soils, and so is the kNN.** Capillary
  flow ignores aggregation, film flow and differences in tortuosity. But
  with the cap, the physics still beats the kNN in every texture group.

---

## Part B. Texture prediction

### B1. What it is for, and the result

The tool reads texture from the reference soils (classifier and neighbour
vote). A method that reads texture from the curve alone, with no reference,
would give an independent second opinion, as the capillary bundle does for
Ks.

**Result: no physical method was good enough.** Every method tested is
significantly worse than the tool on the class. None helps when averaged
with the tool. Some give a weak confidence signal (B2, B3). So the tool reports
no physical texture estimate.

### B2. Four reference-free methods

Four ways to read texture from the retention curve, with no reference data,
were tested as a physical second opinion on the class and the fractions:

- **Inverse Arya–Paris.** The capillary bundle, run toward texture. It
  inverts Arya and Paris (1981): pore radius → particle radius, with their
  scaling exponent α = 1.38, then cumulative particle-size distribution →
  USDA fractions and class. It does no fitting.
- **Inverse Saxton et al. (1986).** These are the retention equations that
  Pan et al. (2010) invert. Pan et al. estimate the mean and spread of
  sand, clay and rock in a field from the spatial mean and spread of
  measured soil moisture. That needs many locations and gives no texture
  for one sample. What carries over to one curve is the forward model: the
  (sand, clay) whose Saxton curve best matches the measured one, on a 1 %
  grid within the stated range of the equations (sand 5–95 %, clay
  5–60 %). Saxton's equations are a regression on measured soils, not
  physics.
- **Fredlund grain-size back-analysis, as read from the paper.** Nguyen and
  Elwood (2025) adjust a Fredlund et al. (2000) unimodal grain-size curve
  and a packing porosity until the Fredlund (2000) physico-empirical
  retention curve matches the measured one. In this reading:
  - each uniform grain fraction gets its own Fredlund and Xing (1994)
    curve;
  - the air entry comes from Sakaki et al. (2014), and the shape from
    Fredlund's effective-diameter correlation;
  - the fractions fill the porosity from the finest up.
- **Fredlund grain-size back-analysis, the authors' code.** This is their
  Python listing, ported with the fewest repairs needed to run it:
  indentation, undefined names, a model that returns the curve at the
  measured suctions, and the sample's own bulk density, porosity and
  saturated water content in place of the values of the listing's test
  soil (Clavet sand). Its algorithm is unchanged:
  - one air entry for the whole soil, from D30;
  - the shape of the curve from the effective diameter of the whole soil;
  - the curve fit adjusts only two grain-size parameters for each random
    draw of the other four;
  - every solution with R² > 0.98 is kept. Here the best one is scored,
    with 100 draws per soil against the listing's 1,500.

  It was run with the listing's silty-sand parameter ranges, and again with
  the paper's fine-grained ranges added.

All four were scored paired against the tool on 1,733 balanced targets
(150 per class, every silt). Each target's whole source was held out from
the tool; the physics needs nothing held out. (2026-10, the 21,298-layer
default reference.)

| | Exact | Group | Macro-F1 | Fraction error (points) |
|---|---|---|---|---|
| Tool (classifier) | 30.8 % | 57.3 % | 28.7 | 13.5 |
| Neighbour vote | 29.2 % | 54.7 % | 28.6 | – |
| Inverse Saxton et al. (1986) | 24.9 % | 54.4 % | 21.6 | 17.9 |
| Fredlund back-analysis, read from the paper | 21.8 % | 51.1 % | 17.9 | 16.0 |
| Inverse Arya–Paris | 21.6 % | 44.6 % | 17.3 | 18.2 |
| Fredlund back-analysis, authors' code (fine ranges added) | 11.7 % | 25.4 % | 6.5 | 31.0 |
| Fredlund back-analysis, authors' code (silty-sand ranges) | 10.3 % | 22.2 % | 3.9 | 30.8 |

Every method is significantly worse than the tool on the class (Saxton
−5.9 pp, p = 4 × 10⁻⁶; the others −9 to −21 pp). Each fails in a different
place:

- **Arya–Paris** almost matches the tool on sandy soils (58 % against
  65 %). It fails on silty and clayey soils (10 % against 25–30 %; clay
  error 15.9 against 8.6 points). There, aggregation and film flow break
  the capillary assumption.
- **Saxton** holds up on clayey soils (28 % against 30 %) but loses on
  sandy ones (45 % against 65 %).
- **The paper-based Fredlund reading** gives the best physical fractions.
  Its clay ranks soils as well as the tool's (Spearman 0.75), and its sand
  slightly better (0.67 against 0.63). It beats the tool on loamy (26 %
  against 18 %) and clayey soils (31 % against 30 %). It fails on sandy
  (20 %) and silty soils (5 %). For 20 % of soils, fits that are almost
  equally good differ by more than 10 points of clay.
- **The authors' code** reads about 80 % of all soils as sand or loamy
  sand. Its median clay is 0 %, against 24 % measured, and clayey soils
  score 0–1 %. The shape of the curve depends only on D30 and on the
  effective diameter of the whole soil, so it carries almost no
  information about the fines. Also, many distributions meet R² (a
  correlation): for a third to a half of the soils, the kept solutions
  differ by more than 20 points of sand. This is the non-uniqueness, and
  the limit to sand and silt, that the authors report.

**They do not help in combination either.** A 50/50 average of any of them
with the tool's fractions raises the error (to 13.6–18.5 points). The class
read from the averaged fractions is lower than the tool's (25.7–28.7 % for
Arya–Paris, Saxton and the paper-based reading).

**As an agreement check, Arya–Paris tells the most.** Where the tool and
its neighbour vote agree, the tool is right 57 % of the time when
Arya–Paris also agrees, and 27 % when it does not (Saxton 47 against 30 %,
the paper-based Fredlund 48 against 33 %). The authors' code splits more
sharply (79 against 31 %, silty-sand ranges). But it agrees with the tool
on only 10 % of the soils, so it rarely says anything. An earlier
Arya–Paris figure of 35.1 % came from the sand-heavy natural class mix of
the reference, and was withdrawn.

### B3. Clay from the dry end

The literature review (2026-10-30) calls this the best-validated inverse
route. At high suction, adsorption on particle surfaces holds the water,
so the water at 1500 kPa scales with clay (Tuller and Or 2005; Wuddivira et
al. 2012). Three forms were tested (`verify_literature_methods.py`):

- the NRCS laboratory rule, clay ≈ 2.5 × gravimetric water at 1500 kPa
  (%), with the water taken from the fitted curve and divided by bulk
  density;
- the same rule, with organic carbon (%) taken off the water first;
- a straight line of clay on that water, or on the volumetric water,
  fitted on the other studies.

They were scored against the tool's fractions on the 1,731 targets of
`verify_predictive.py`. There are 71 studies; each study is hidden from the
tool, and depth and bulk density are supplied. The tool's predictions are
read from `figures/cache/predictive.csv`.

| Clay estimate | n | Mean error (points) | Bias | ρ | vs tool | Studies better |
|---|---|---|---|---|---|---|
| **Tool (kNN fractions)** | 1,730 | **8.45** | −0.92 | 0.76 | – | – |
| NRCS rule | 1,730 | 9.95 | +4.98 | 0.77 | +1.50 (p < 10⁻¹¹) | 33/71 |
| NRCS rule, organic carbon removed | 741 | 8.62 | +1.50 | 0.75 | +0.36 (p = 0.60) | 27/54 |
| Line on volumetric water, other studies | 1,730 | 9.01 | −1.99 | 0.76 | +0.56 (p < 10⁻⁸) | 35/71 |
| Line on gravimetric water, other studies | 1,730 | 10.42 | −3.48 | 0.74 | +1.97 (p < 10⁻²⁹) | 22/71 |
| Half tool, half gravimetric line | 1,730 | 8.95 | −2.20 | 0.77 | +0.50 (p < 10⁻⁷) | 27/71 |

- **No form beats the tool.** Only the rule with organic carbon removed is
  level with it, and only 741 of the targets have organic carbon.
- **By measured texture group, the rule beats the tool in clayey soils
  only** (11.7 against 12.3 points, p = 0.03). It loses in sandy (6.9
  against 3.9), loamy (9.6 against 7.8) and silty soils (10.8 against
  8.6). A user does not know the group in advance, so this cannot be used
  directly.
- **The class gets worse.** Keep the tool's sand and replace its clay with
  the fitted line: the class read from the fractions falls from 29.5 to
  23.6 % exact (27.0 % with half and half). The tool's own classifier
  scores 32.5 % on these targets.
- **As a confidence flag for the fractions, it has some value.** Where the
  rule (organic carbon removed when known) agrees with the tool's clay
  within 5 points (54 % of soils), the tool's clay error is 7.1 points.
  Where it does not agree, the error is 10.0. It does not flag the class
  (32.6 against 32.3 % exact). Not adopted.
- **The test is partly circular.** Where GSHP had no dry-end points, its θr
  came from a pedotransfer function on texture (Gupta et al. 2022). That
  affects the rule and the tool's classifier equally, because both read
  the same fitted curve at 1500 kPa.

### B4. Not tested, and why

The other texture routes in the review were not tested, for two reasons:

- **They belong to a family already scored here.** Mohammadi and
  Vanclooster, Arya and Heitman, and the fractal and lognormal models
  refine inverse Arya–Paris, and they have the same documented failure in
  fine soils.
- **They need data that one curve does not have:** humidity isotherms
  below −10 MPa, or time series of soil moisture for data assimilation.

Also, any method that reads the same curve is bounded by the Cover and
Hart ceiling: about 44 % exact for a new source (README, *Validation*).

---

## Reproduce

Run from the repository root. Set one thread per numerical library (the
runs are single-threaded by design). Use `caffeinate -i` so that macOS
does not sleep during the runs.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
caffeinate -i .venv/bin/python verify_ks_physical.py > ks_physical.log 2>&1              # ~10 min: the source-blocked tables of A3
caffeinate -i .venv/bin/python verify_external.py nj_ssir26 > nj_external.log 2>&1        # the New Jersey test set, all three designs
caffeinate -i .venv/bin/python figures/make_external_report.py nj_ssir26                  # its fig6 and ks_variants.csv
caffeinate -i .venv/bin/python figures/make_external_report.py babaeian_zanjanrood        # the same for Zanjanrood
.venv/bin/python verify_literature_methods.py > literature_methods.log 2>&1               # seconds: the literature forms (A4 and B3); reads figures/cache/ks_neighbors.csv and predictive.csv
```

The per-soil reports go to `figures/external/<set>/`. `caffeinate` is for
macOS only; on other systems, leave it out.

## References

- Ahuja, L.R., Cassel, D.K., Bruce, R.R. and Barnes, B.B. (1989).
  Evaluation of spatial distribution of hydraulic conductivity using
  effective porosity data. *Soil Science* 148:404–411.
  doi:10.1097/00010694-198912000-00002.
- Arya, L.M. and Paris, J.F. (1981). A physicoempirical model to predict the
  soil moisture characteristic from particle-size distribution and bulk
  density data. *Soil Science Society of America Journal* 45:1023–1030.
- Bouwer, H. (1966). Rapid field measurement of air entry value and
  hydraulic conductivity of soil as significant parameters in flow system
  analysis. *Water Resources Research* 2:729–738.
- Brakensiek, D.L., Engleman, R.L. and Rawls, W.J. (1981). Variation within
  texture classes of soil water parameters. *Transactions of the ASAE*
  24:335–339.
- Brooks, R.H. and Corey, A.T. (1964). *Hydraulic properties of porous
  media*. Hydrology Paper 3, Colorado State University, Fort Collins.
- Brutsaert, W. (1967). Some methods of calculating unsaturated
  permeability. *Transactions of the ASAE* 10:400–404.
- Childs, E.C. and Collis-George, N. (1950). The permeability of porous
  materials. *Proceedings of the Royal Society of London A* 201:392–405.
- Fredlund, D.G. and Xing, A. (1994). Equations for the soil-water
  characteristic curve. *Canadian Geotechnical Journal* 31:521–532.
- Fredlund, M.D. (2000). *The role of unsaturated soil property functions in
  the practice of unsaturated soil mechanics*. PhD thesis, University of
  Saskatchewan, Saskatoon.
- Fredlund, M.D., Fredlund, D.G. and Wilson, G.W. (2000). An equation to
  represent grain-size distribution. *Canadian Geotechnical Journal*
  37:817–827.
- Ghanbarian, B., Hunt, A.G., Ewing, R.P. and Skinner, T.E. (2014).
  Universal scaling of the formation factor in porous media derived by
  combining percolation and effective medium theories. *Geophysical
  Research Letters* 41:3884–3890.
- Ghanbarian, B. and Patrignani, A. (2026). Estimating saturated hydraulic
  conductivity of soil samples from Kansas Mesonet database. *Soil Science
  Society of America Journal*. doi:10.1002/saj2.70222. Preprint:
  arXiv:2602.00893.
- Greenland, D.J. (1977). Soil damage by intensive arable cultivation:
  temporary or permanent? *Philosophical Transactions of the Royal Society
  of London B* 281:193–208.
- Guarracino, L. (2007). Estimation of saturated hydraulic conductivity Ks
  from the van Genuchten shape parameter α. *Water Resources Research*
  43:W11502. doi:10.1029/2006WR005766.
- Gupta, S., Papritz, A., Lehmann, P., Hengl, T., Bonetti, S. and Or, D.
  (2022). Global Soil Hydraulic Properties dataset based on legacy site
  observations and robust parameterization. *Scientific Data* 9:444.
- Hillel, D. (1980). *Fundamentals of Soil Physics*. Academic Press, New
  York. Chapters 8 and 9.
- Ippisch, O., Vogel, H.-J. and Bastian, P. (2006). Validity limits for the
  van Genuchten–Mualem model and implications for parameter estimation and
  numerical simulation. *Advances in Water Resources* 29:1780–1789.
- Jackson, R.D. (1972). On the calculation of hydraulic conductivity. *Soil
  Science Society of America Proceedings* 36:380–382.
- Johnson, D.L., Koplik, J. and Schwartz, L.M. (1986). New pore-size
  parameter characterizing transport in porous media. *Physical Review
  Letters* 57:2564–2567.
- Marshall, T.J. (1958). A relation between permeability and size
  distribution of pores. *Journal of Soil Science* 9:1–8.
- Mishra, S. and Parker, J.C. (1990). On the relation between saturated
  conductivity and capillary retention characteristics. *Groundwater*
  28:775–777.
- Nguyen, M. and Elwood, D. (2025). An iterative program to back-analyze
  grain-size distribution from a predetermined soil–water characteristic
  curve. *Transportation Geotechnics* 51:101485.
- Pan, F., Peters-Lidard, C.D. and King, A.W. (2010). Inverse method for
  estimating the spatial variability of soil particle size distribution
  from observed soil moisture. *Journal of Hydrologic Engineering*
  15:931–938.
- Peters, A., Hohenbrink, T.L., Iden, S.C., van Genuchten, M.Th. and
  Durner, W. (2023). Prediction of the absolute hydraulic conductivity
  function from soil water retention data. *Hydrology and Earth System
  Sciences* 27:1565–1582.
- Sakaki, T., Komatsu, M. and Takahashi, M. (2014). Rules-of-thumb for
  predicting air-entry value of disturbed sands from particle size. *Soil
  Science Society of America Journal* 78:454–464.
- Saxton, K.E., Rawls, W.J., Romberger, J.S. and Papendick, R.I. (1986).
  Estimating generalized soil-water characteristics from texture. *Soil
  Science Society of America Journal* 50:1031–1036.
- Skaggs, T.H. (2011). Assessment of critical path analyses of the
  relationship between permeability and electrical conductivity of pore
  networks. *Advances in Water Resources* 34:1335–1342.
- Tuller, M. and Or, D. (2005). Water films and scaling of soil
  characteristic curves at low water contents. *Water Resources Research*
  41:W09403. doi:10.1029/2005WR004142.
- van Genuchten, M.Th. (1980). A closed-form equation for predicting the
  hydraulic conductivity of unsaturated soils. *Soil Science Society of
  America Journal* 44:892–898.
- Vogel, T., van Genuchten, M.Th. and Císlerová, M. (2001). Effect of the
  shape of the soil hydraulic functions near saturation on variably-saturated
  flow predictions. *Advances in Water Resources* 24:133–144.
- Wuddivira, M.N., Robinson, D.A., Lebron, I. et al. (2012). Estimation of
  soil clay content from hygroscopic water content measurements. *Soil
  Science Society of America Journal* 76:1529–1535.
  doi:10.2136/sssaj2012.0034.
