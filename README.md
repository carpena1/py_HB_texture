# py_HB_texture — soil texture and Ks from a measured water retention curve

Give it a measured soil water characteristic curve (SWCC) — from an in-situ
sensor or a laboratory core — and it returns the USDA texture class, the
sand/silt/clay fractions and the saturated hydraulic conductivity Ks, each
with an uncertainty range and the real reference soils it was matched to.

The tool fits a van Genuchten curve to the data and compares it with about
21,000 measured soil layers from nine databases. It is built for curves
measured on undisturbed soil — typically in-situ sensors — and the reference
can grow with your own lab-verified sites.

## What to expect

The tool always searches the whole reference. How well it does depends on
whether curves from your data source — your laboratory, or your sensor
network — are already in it, because curves from one source share a
measurement signature about as strong as the texture signal:

| | new data source | your source already in the reference | chance |
|---|---|---|---|
| exact USDA class (12 classes) | **29 %** (best possible* 42 %) | **37 %** (best possible 53 %) | 8 % |
| … with depth and bulk density | 29 % | 43 % | |
| texture group (4 groups) | **56 %** (68 %) | **63 %** (73 %) | 25 % |
| … with depth and bulk density | 57 % | 67 % | |
| Ks, typical error | factor of 8 | factor of 3 | — |
| … with bulk density | factor of 7 | factor of 2.6 | — |
| Ks within a factor of 10 | 53 % | 76 % | — |

\* The Cover & Hart (1967) bound for *any* method that uses the four van
Genuchten parameters (see "Validation"). All rows are the same 1,733 test
soils (`verify_holdout.py`); the "with" rows use the 1,681 of them that carry
depth and bulk density, compared with the curve alone on those soils.

A new depth at a site whose other horizons are already in the reference does
a little better than a new site (39.5 % against 37.2 %, p=0.012). Depth and
bulk density add about 5 points of exact class at both (p<0.001), and
nothing for a new source.

**You can move from the first column to the second.** Adding lab-verified
curves from your own sites (see "Adding your own verified data") is the most
effective improvement measured: with a source's own soils in the reference,
texture gains 5 points and Ks ranking goes from nearly useless (ρ = 0.11) to
modest (ρ = 0.36). With those in, supplying depth and bulk density together
adds about 4 points more.

A retention curve does not pin texture down: soils of neighbouring classes
produce nearly identical curves, and the laboratory that measured a curve
leaves a mark on it that is as large as the texture signal. The tool therefore
reports ranked class probabilities, fraction ranges and a Ks band rather than
a single answer, plus a second opinion whose agreement is a useful confidence
flag (right 36 % of the time when the two agree, 24 % when they do not, for a
new source). Treat Ks as an order-of-magnitude estimate.

## Setup

    python3 -m venv .venv
    .venv/bin/pip install -r requirements.txt

Calling `.venv/bin/python` directly (as below) runs inside the virtual
environment without activating it.

The reference tables are committed, so **the tool runs out of the box with no
downloads.** Two tables are exceptions: EU-HYDI, which is restricted to its
consortium, and the Laikipia soils, which are unpublished. Their tables
(`data/euhydi_reference.csv`, `data/willard_reference.csv`) are built locally
with `prepare_euhydi.py` and `prepare_willard.py` by those who hold the data
and are never distributed. A fresh clone runs on the other seven databases
(14,372 layers) and prints a note saying so.

## Usage

    .venv/bin/python swcc_texture.py mydata.csv

`mydata.csv` has two columns (comma or whitespace separated): suction h in
kPa and volumetric water content θ, as a fraction or in % (detected
automatically). At least five points; a curve that reaches both near
saturation and the dry end (≥ 1,000 kPa) carries the most information.

Options a user needs:

    --sample-type undisturbed   how the curve was measured (default). Use it
                                for in-situ sensors and intact cores.
    --sample-type disturbed     repacked or sieved samples
    --sample-type unknown       no information
    --depth 30                  sample mid-depth in cm (optional)
    --bulk-density 1.35         oven-dry bulk density in g/cm3 (optional)
    --andic yes                 the soil has andic properties (optional;
                                --andic no when known not to)
    --json out.json             also write the full result as JSON

Supplied together, depth and bulk density add about 4 points of exact class
when your own verified data are in the reference; for a data source the
reference has never seen their effect is small and not stable. Supply them
when you have them, knowing that part of the gain favours the classes your
own data contain most (see "Local covariates"). Bulk density also joins the
neighbour search, which lowers the Ks error slightly.

`--andic yes` is for soils known to have andic properties — an Andisol or
Andosol, or andic by oxalate Al + ½Fe and phosphate retention — typically on
volcanic ash. For such soils from a source the reference has never seen it
adds 9 to 14 points of exact class and 12 to 21 of texture group, and brings
the fractions 2 points closer; Ks does not change (see "Volcanic and andic
soils"). Leave it out when you do not know.

From Python: `swcc_texture.estimate(h, theta, sample_type="undisturbed")`;
for depth and bulk density, pass a classifier trained with them
(`TextureGBM(covariates=["depth_cm", "bd"])`) as `clf=`, a reference built
with `GshpReference(use_bd=True)` as `ref=`, and `depth=` and
`bulk_density=`. For andic properties, pass `TextureGBM(covariates=["andic_code"])`
and `andic="yes"` or `"no"`.
`--help` lists further options; they exist for the verification scripts, and
the defaults are the combination that tested best.

### Reading the output — three worked examples

**1. Silt loam — a real soil, read well** (`Testing/13_SiL_UNSODA.csv`).
The measured curve of a Humbeek silt loam (Belgium, Ap horizon, 0–30 cm,
intact core at the wet end; UNSODA code 4070), a soil that is not in the
default reference (it is in `data/unsoda_reference.csv`, which only
`--reference all` loads):

    $ .venv/bin/python swcc_texture.py Testing/13_SiL_UNSODA.csv

    van Genuchten fit (m = 1-1/n, alpha in kPa^-1):
      thetar = 0.0672  thetas = 0.3990  alpha = 0.0400  n = 1.649  m = 0.394  (RMSE 0.0074)

    Predicted USDA texture class: SILT LOAM   (from the gradient-boosted classifier)
      silt loam         51.6 %
      loam              16.0 %
      sandy loam         9.7 %
      sandy clay loam    6.0 %
      loamy sand         3.8 %
      kNN second opinion agrees: silt loam (56.5 %)

    Particle fractions, mean [5-95 %]  (from the kNN; mean plots as: silt loam):
      sand   32.6 %  [  8.7 -  66.1]
      silt   50.4 %  [ 17.7 -  74.6]
      clay   17.0 %  [  8.5 -  29.0]

    Ks = 2.15 cm/h  [5-95 %: 0.0171 - 36.2]  (from ~15 of 30 undisturbed neighbors with measured Ksat)

- **A clear answer.** The classifier puts more than half its probability on
  the true class, three times the next one, and the neighbour vote agrees —
  the combination to trust most, though not blindly: when the two agree the
  class is right about half the time (see "What to expect").
- **Fractions close to the laboratory's** (35 / 55 / 10 % sand / silt /
  clay): within 3 points for sand, 5 for silt and 7 for clay, all inside
  their ranges.
- **Ks within a factor of two** (2.15 against 1.21 cm/h measured), taken
  from undisturbed neighbours only because the wet end was measured on an
  intact core.

The next two curves are generated from the Carsel & Parrish (1988)
class-typical parameters: class averages rather than real samples, so they
are run with `--sample-type unknown` (real sensor data should keep the
default). They show what the tool does when a curve is ambiguous, and when it
is confidently wrong.

**2. Sandy clay loam — a disagreement that points to the right answer**
(`Testing/7_SCL_CnP.csv`):

    $ .venv/bin/python swcc_texture.py Testing/7_SCL_CnP.csv --sample-type unknown

    van Genuchten fit (m = 1-1/n, alpha in kPa^-1):
      thetar = 0.1000  thetas = 0.3900  alpha = 0.5900  n = 1.480  m = 0.324  (RMSE 0.0000)

    Predicted USDA texture class: SANDY LOAM   (from the gradient-boosted classifier)
      sandy loam        25.7 %
      sandy clay loam   20.4 %
      loam              13.1 %
      loamy sand        11.1 %
      sandy clay         8.6 %
      kNN second opinion DISAGREES: sandy clay loam (43.9 %)

    Particle fractions, mean [5-95 %]  (from the kNN; mean plots as: sandy loam):
      sand   68.7 %  [ 12.0 -  84.0]
      silt   14.6 %  [  5.0 -  67.0]
      clay   16.7 %  [  5.3 -  25.0]

    Ks = 4.71 cm/h  [5-95 %: 0.00642 - 21.2]  (from ~17 of 30 neighbors with measured Ksat)

- **A weak preference.** The classifier gives its top class only 25.7 % and
  spreads the rest over neighbouring classes; the true class, sandy clay
  loam, is second (20.4 %).
- **The second opinion disagrees, and it is right.** The neighbour vote names
  sandy clay loam with 43.9 %. When the two disagree, the classifier's answer
  is right only about one time in five for a new source: read a disagreement
  as "one of these two, check both".
- **The fraction mean plots as sandy loam.** The mean of a set of neighbour
  fractions need not fall in the modal class; the intervals are the better
  guide.
- **Ks is inside its band but 3.6 times too high** (4.71 against 1.31 cm/h),
  and the band spans more than three orders of magnitude.

**3. Loam — a confident wrong answer that the agreement flag does not catch**
(`Testing/4_L_CnP.csv`):

    $ .venv/bin/python swcc_texture.py Testing/4_L_CnP.csv --sample-type unknown

    Predicted USDA texture class: LOAMY SAND   (from the gradient-boosted classifier)
      loamy sand        49.3 %
      loam              13.6 %
      sandy loam        10.3 %
      silt loam          6.0 %
      sandy clay loam    5.7 %
      kNN second opinion agrees: loamy sand (49.0 %)

    Particle fractions, mean [5-95 %]  (from the kNN; mean plots as: sandy loam):
      sand   73.2 %  [ 19.3 -  89.7]
      silt   19.3 %  [  2.0 -  68.6]
      clay    7.5 %  [  1.0 -  22.1]

    Ks = 7.65 cm/h  [5-95 %: 2.1 - 25.9]  (from ~23 of 30 neighbors with measured Ksat)

- **Wrong, confidently, and both halves agree on it.** Half the probability
  goes to loamy sand, and the neighbour vote says the same with the same
  weight; the true class, loam, is second (13.6 %). Agreement raises the odds
  of being right but does not guarantee it: the Carsel & Parrish loam sits
  where the reference holds loamy sands (the class-mean benchmarks below
  explain why).
- **The warnings are in the fractions.** Their mean (73 / 19 / 8 % sand /
  silt / clay, against the class centroid's 41 / 40 / 18) plots as sandy loam,
  not as the predicted class, and the silt range runs from 2 to 69 %.
- **Ks is 7 times too high, and the measured value (1.04 cm/h) falls below the
  whole band** — a wrong texture neighbourhood gives a wrong Ks.

## Adding your own verified data

When some of your sites are sent to the laboratory for texture (and, ideally,
Ks measured on undisturbed cores), feed them back into the reference. Later
predictions from the same sensors or laboratory are then matched against your
own soils — the "already in the reference" column above.

Prepare one CSV with one row per retention point. Every row needs only
`sample_id`, `h_kpa` and `theta`; the per-sample values (texture, Ks, depth,
coordinates, ...) need to be written once, on any row of the sample, and the
other rows can stop after `theta` or leave those fields empty. Repeating them
on every row also works; if they differ, the first non-empty value is used.

    sample_id,h_kpa,theta,sand,silt,clay,ksat_cmh,depth_cm,bd,sample_type,site,lat,lon
    S01-30,0,0.43,60,13,27,1.3,30,1.40,undisturbed,S01,29.64,-82.35
    S01-30,5,0.36
    S01-30,33,0.21
    ...
    S02-30,0,0.47,20,55,25,,30,1.30,undisturbed,S02,29.71,-82.40
    S02-30,5,0.44
    ...

`sample_id`, `h_kpa`, `theta`, `sand`, `silt` and `clay` are required;
`ksat_cmh` (cm/h), `depth_cm`, `bd` (g/cm3), `sample_type` (default
undisturbed), `site` (groups the samples of one location) and `lat`/`lon` (decimal
degrees, WGS84; south and west negative) are optional. Then:

    .venv/bin/python add_local_data.py mycurves.csv --source mynetwork

Each sample is fitted with the tool's own routine (at least five points, fit
RMSE ≤ 0.03) and written to `data/local_reference.csv`, which the default
reference picks up automatically; adding a `sample_id` again replaces it. The
file is git-ignored, so your data stays private. The more of your texture
range the verified sites cover, the better — a few sites per soil type are
worth more than many of one.

## How it works

1. **Fit.** Van Genuchten parameters (θr, θs, α, n) are fitted by least squares
   with the Mualem constraint m = 1 − 1/n. The fit's covariance is sampled
   300 times, so fit uncertainty flows into every output.
2. **Texture class — gradient-boosted classifier.** Trained on the reference
   layers' four van Genuchten parameters plus sample type (and depth and bulk
   density when you supply them), class-balanced; its probabilities are
   averaged over the 300 draws.
3. **Everything else — nearest neighbours.** For each draw the 30 nearest
   reference layers in standardised parameter space (plus bulk density, when
   you supply it) vote, weighted by 1/d² and
   by inverse class frequency. They supply the particle fractions and their
   ranges, the list of matched soils, and the second-opinion class. **Ks** comes
   from the neighbours that carry a measured Ks; for an undisturbed sample,
   only from *undisturbed* neighbours, because near saturation an intact core
   keeps the large pores that repacking destroys. Disturbed and unknown
   samples use all soils, since the reference holds only 45 disturbed layers
   with a measured Ks.

The classifier and the neighbours see the same draws. Only the class comes
from the classifier; `--model knn` swaps it for the neighbour vote and leaves
every other number identical.

## Reference data

| source | layers | with Ks | sample type | distribution |
|---|---|---|---|---|
| GSHP (Gupta et al. 2022), quality-filtered | 9,996 | 6,644 | per contributor | CC BY 4.0, committed |
| NCSS/KSSL (USDA-NRCS) | 2,530 | 0 | undisturbed | public domain, committed |
| Hohenbrink et al. (2023), German cores | 560 | 402 | undisturbed | committed |
| Zanjanrood, north-west Iran (Babaeian et al. 2015) | 169 | 169 | undisturbed | committed |
| Arizona, USA (E. Babaeian, measurements shared with this project) | 21 | 21 | undisturbed | committed |
| Canary Islands andic soils (Armas Espinel 2013) | 66 | 62 | undisturbed | committed |
| Yellow River Basin, China (Tong et al. 2024) | 1,030 | 1,015 | undisturbed | committed |
| EU-HYDI (Weynants et al. 2013), 21 laboratories, 16 countries | 6,797 | 2,704 | per sample | **restricted, never distributed** |
| Laikipia, central Kenya (Willard et al., under review) | 86 | 0 | undisturbed | **unpublished, never distributed** |
| **default reference** | **21,255** | **11,017** | | |

Each source gets its own `prepare_*.py`, because each ships something
different, and each needed a correction before its curves were comparable:

- **KSSL** has no reading wetter than 6 kPa, so nothing constrains saturation;
  a plain fit puts α five to sixteen times too low. One anchor point at
  θ(0) = 0.95 × porosity fixes it (θs sits ~5 % below porosity because of
  entrapped air; median θs/porosity in GSHP is 0.946). KSSL has no Ks.
- **Hohenbrink** has ~199 HYPROP points on the wet limb against three WP4
  points on the dry one; points are thinned to one median per pF bin so the
  fit is not weighted 98:2 away from the part of the curve that carries
  texture.
- **EU-HYDI** samples need at least five reliable retention points and USDA
  fractions; Ks is taken only from saturated-conductivity methods, and curves
  are thinned in pF as for Hohenbrink. Several contributors reported only two
  or four points per sample (Greece, most of Spain), too few to fit a curve.
  Polish fractions are summed from raw size classes breaking at 2 and 50 µm.
- **Zanjanrood** curves combine an intact core for 0–100 cm suction with
  sieved soil for the dry end, and are refitted here; four curves with no
  point wetter than 330 cm are dropped. Ks was measured on repacked samples
  but is used as undisturbed: predicted from the rest of the reference before
  the merge, it came out at the right level (`verify_external.py`).
- **Laikipia** clays (74 of 86 samples) on the northern slopes of Mount Kenya
  were measured on a sandbox to 10 kPa and pressure plates beyond. The
  sandbox did not saturate the cores, and the correction applied to it (a
  constant raise of pF 0–2 to 95 % of porosity) left a step between 10 and
  14 kPa holding half of each curve's water loss; read as delivered, the tool
  took the clays for silty clay loam and loamy sand. Only saturation, at 95 %
  of porosity, and the pressure-plate points are used. Ks was measured only in
  the field, in the dry season, when these vertic clays were cracked (median
  178 cm/h on the drier clay-rich samples), so none of it enters the
  reference; the field values stay in the table.
- **Arizona** curves run from saturation to about 29,000 kPa; only points to
  1500 kPa are fitted, as for the rest of the reference, and a 500 cm reading
  that falls below the 800 cm one in three samples (AZ15, AZ18, AZ19) is
  replaced by interpolation between its neighbours. There are no coordinates
  or depths, so the soils are not on the map.
- **Canary Islands** retention points are not tabulated in the thesis; they
  are read from its vector figures (Tempe cells to 90 kPa, pressure plates at
  100, 500 and 1500 kPa) and check against its tables (water at 1500 kPa
  within 0.3 points) and its own van Genuchten fits. Texture is the
  hydrometer analysis after hexametaphosphate dispersion, which under-disperses
  andic material; the resin-dispersed fractions are kept alongside.
- **Yellow River Basin** curves (centrifuge, 1–1000 kPa, 475 sites to 5 m
  deep) are published only as the authors' van Genuchten fits, and in 170 of
  them the fitted θs is more than 0.05 from the measured one. Each is rebuilt
  at ten suctions from its fit, the measured θs is added at h = 0, and the
  whole is refitted here. Texture is by laser diffraction, which reads about
  half the clay of the pipette and hydrometer methods elsewhere in the
  reference (median 10 % against 18 % for silt loam). As a new source the set
  reads poorly (14 % exact class, `verify_external.py`): with its large α the
  curves look coarser than they are. Added to the reference it does no harm
  to the other soils (+1.3 pp exact, p=0.15) and helps silt loam (+6 pp,
  p=0.08).

**Sample type** (`sample_type`, with its evidence in `sample_type_source`) is
set by how the wet end of each curve was measured: an intact core for the wet
end with sieved material for the dry end — standard practice — counts as
undisturbed; disturbed means the wet end itself was measured on repacked soil.
Labels come from each source's metadata: GSHP's own field, EU-HYDI's method
table per sample, the KSSL method codes and the Hohenbrink data description.
Florida's soils (5,735 layers), which GSHP lists as unknown, are undisturbed:
one of the dataset's authors recalls the rings for conductivity and water
release being taken in the field for each horizon, without repacking. In the
default reference: 18,453 undisturbed, 1,581 disturbed, 1,221 unknown. Of
the 11,017 layers with a measured Ks, 10,951 are undisturbed and only 45
disturbed.

**Coverage.** 20,598 layers carry coordinates (`figures/fig2_coverage_map.png`).
North America holds 9,256, Europe 7,740 (38 %), Asia 2,653 (1,030 of them in
the Yellow River Basin), South America 639, Africa 750 (with the Canary
Islands) and Australasia 111; Florida alone supplies 28 %. Southern Europe and the Balkans remain thin.

**Texture classes.** The classes are far from balanced
(`figures/fig5_texture_classes.png`): sand outnumbers silt 56 to 1, and
Florida alone supplies 79 % of the sand (3,657 of 4,620 layers). The
restricted tables supply most of the silty clay (57 %) and 40 % of the silt
loam; the Yellow River Basin adds 732 silt loams and 22 of the 83 silts.
Ks is measured for 80 % of sandy soils but for about a third of silty ones.

| group | class | layers | share | with Ks | distributed | restricted |
|---|---|---|---|---|---|---|
| Sandy | sand | 4,620 | 21.7 % | 4,020 | 4,051 | 569 |
| | loamy sand | 1,467 | 6.9 % | 856 | 1,000 | 467 |
| Loamy | sandy loam | 3,788 | 17.8 % | 1,726 | 2,306 | 1,482 |
| | loam | 1,882 | 8.9 % | 516 | 1,162 | 720 |
| | sandy clay loam | 1,323 | 6.2 % | 779 | 1,101 | 222 |
| | clay loam | 1,007 | 4.7 % | 370 | 632 | 375 |
| Silty | silt | 83 | 0.4 % | 48 | 69 | 14 |
| | silt loam | 3,328 | 15.7 % | 1,331 | 1,989 | 1,339 |
| | silty clay loam | 1,326 | 6.2 % | 327 | 721 | 605 |
| Clayey | sandy clay | 185 | 0.9 % | 112 | 172 | 13 |
| | silty clay | 627 | 2.9 % | 182 | 270 | 357 |
| | clay | 1,619 | 7.6 % | 750 | 899 | 720 |

| texture group | layers | share | with Ks | distributed | restricted |
|---|---|---|---|---|---|
| Sandy | 6,087 | 28.6 % | 4,876 (80 %) | 5,051 | 1,036 |
| Loamy | 8,000 | 37.6 % | 3,391 (42 %) | 5,201 | 2,799 |
| Silty | 4,737 | 22.3 % | 1,706 (36 %) | 2,779 | 1,958 |
| Clayey | 2,431 | 11.4 % | 1,044 (43 %) | 1,341 | 1,090 |
| **total** | **21,255** | | **11,017** | **14,372** | **6,883** |

The neighbour vote weights each class by the inverse of its count and the
classifier uses balanced class weights, so neither favours sand for being
common. Weighting cannot add information, though: silt (83 layers) and sandy
clay (185) remain the least-represented classes, and the per-class accuracy
in "Validation" reflects it.

UNSODA 2.0 (581 soils) and the Belgian sDB (165) are also prepared and can be
loaded with `--reference all`, but most of their soils are already inside
GSHP, so they are not in the default.

**Volcanic parent material and andic properties** (`volcanic`, `andic`,
joined from `data/volcanic_flags.csv` by `prepare_volcanic.py`) are two
separate fields, because they differ: the Laikipia soils sit on Mount
Kenya's volcanics but show no andic properties. Evidence is taken from the
sources first — KSSL taxonomy and andic lab criteria, the Africa Soil
Profiles database's parent material and class, the EU-HYDI report's own
description of each contributor's sites — then from the SoilGrids 2.0
Andosols probability and the distance to Smithsonian GVP volcanoes, which
are specific but insensitive (checked against the KSSL labels). In the
default reference: 524 layers are known volcanic and 42 probable; 252 are
andic and 155 likely andic (Campania, the Canary Islands; Kamchatka, from
bulk density ≤ 0.90).

## Validation

Accuracy is measured for both kinds of user by hiding known soils and
predicting them. For a **new source**, each target's whole contributing
laboratory is removed from the reference during the test; for a source
**already in the reference**, the target's own profile, or only the target
layer itself. In real use nothing is removed — this is only how the
situations are simulated.
Matching a soil against other soils from its own laboratory is worth about 10
points of exact class. All figures come from paired tests on the same
targets, with the classifier and the neighbours both blind to the target.

### Texture

The same 1,733 targets (150 per class, every silt) under each hold-out, with
the shipped tool (`verify_holdout.py`, `figures/validation/v6_holdout_levels.png`).
The last two columns run it as with `--depth` and `--bulk-density`, on the
1,681 targets that carry both:

| hidden along with the test soil | real-life situation | exact class | texture group | exact, + depth & BD | group, + depth & BD |
|---|---|---|---|---|---|
| nothing else | a new depth at a site already in the reference | 39.5 % | 63.6 % | 44.6 % (+5.2, p<0.001) | 66.8 % (+3.0, p=0.004) |
| its profile | a new site, from a source already in the reference | 37.2 % | 62.7 % | 42.8 % (+5.6, p<0.001) | 66.9 % (+3.9, p<0.001) |
| its whole source | a new site, from a new source | **28.7 %** | **56.3 %** | 29.0 % (+0.2, p=0.85) | 57.0 % (+0.7, p=0.55) |
| best possible (ceiling) | | 53 % / 42 % | 73 % / 68 % | | |

(The gains in brackets are paired against the curve alone on the same 1,681
soils, where it scores 39.4, 37.2 and 28.7 % exact.)

76 of the targets are Yellow River Basin soils, including 22 of the 83 silts,
and they read poorly (see "Reference data"). On the previous target draw,
without them, the current reference scores 40.9, 38.5 and 31.9 % exact —
within a point of the reference before the Yellow River soils joined (41.4,
39.7 and 30.9 %) — so the lower figures above reflect harder targets, not a
worse reference. Part of the depth-and-bulk-density gain with the source in
the reference is the Yellow River soils' own signature (deep loess cores); on
the previous draw it is +3.1 pp for a new site (p=0.003).

Keeping the target's other horizons in the reference helps a little
(+2.3 pp, p=0.012). 1,123 of the targets have a horizon of the same class in
their profile (43.7 % with it in the reference against 40.7 % without);
siblings of a different class neither help nor mislead (32.5 % against
31.3 %) — eluviation makes the horizons of one profile genuinely different
samples. What matters is the data source: hiding it costs 8.5 points
(p<0.0001).

**The ceiling.** Cover & Hart (1967) bound the Bayes error by the
nearest-neighbour error, giving a model-free limit for any classifier on these
four parameters (`verify_ceiling.py`): 1NN accuracy 30.3 % gives a bracket of
[30.3 %, 53.2 %] in-distribution and [20.9 %, 42.3 %] with the laboratory held
out. It is a property of the inputs, not of the database — more reference
data does not raise it, only more informative inputs can.

**Classifier vs neighbour vote.** The gradient-boosted classifier is within
1.5 points of the neighbour vote at both levels (38.6 % against 37.3 % for a
new site, p=0.26; 28.7 % against 27.3 % for a new source, p=0.22); its edge
is the texture group for a new source (56.3 % against 54.6 %,
`verify_hybrid.py`; on the previous reference it was 3 points ahead, p=0.004),
and the neighbour vote supplies the second opinion. Its
weak class is silt, which it rarely predicts.

**Benchmarks** — synthetic curves from class-typical parameters, and real
GSHP soils held out one at a time (`verify_carsel_parrish.py`,
`verify_rosetta.py`, `verify_gshp.py`, `verify_groups.py`):

| benchmark | exact class | true class in top 2 | texture group | Ks inside 5-95 % band | typical Ks error |
|---|---|---|---|---|---|
| Carsel & Parrish (1988) class means | 1/12 | 6/12 | 6/12 | 9/12 | ×4.4 |
| ROSETTA (Schaap et al. 2001) class means | 3/12 | 5/12 | 7/12 | 9/12 | ×4.1 |
| GSHP soils, leave-one-out | 4/12 | 7/12 | 7/12 | 11/12 | ×2.6 |

Twelve soils per benchmark is small; a one-soil change is noise. The class-mean
benchmarks score low partly because the databases disagree on where each class
sits in van Genuchten space (clay α is 0.08, 0.15 and 0.41 kPa⁻¹ in Carsel,
ROSETTA and GSHP), and the fine classes — especially the Clayey group — are
the persistent failure. The GSHP test is in-distribution: each soil's own
laboratory stays in the reference.

### Ks

Ks is harder than texture, because near saturation it depends on structure —
macropores and aggregation — that a texture-matched lookup cannot see. Ks
error for undisturbed soils (leave-one-soil-out: 585 soils,
`figures/make_validation_figures.py`; the other two rows: 480 soils from 21
laboratories, `verify_sample_type.py`):

| | median error | within 5× / 10× | 5–95 % band contains it | rank correlation ρ |
|---|---|---|---|---|
| new depth, other horizons in the reference | 0.39 dex (×2.5) | 68 % / 78 % | 79 % | 0.70 |
| new site, source in the reference | 0.56 dex (×3.6) | 58 % / 74 % | 79 % | 0.64 |
| new site, new source | 0.87 dex (×7.4) | 42 % / 55 % | 68 % | 0.31 |

The band is nominally 90 %; treat it as a minimum range. Supplying bulk
density (`--bulk-density`) tightens the first row to 0.36 dex (×2.3), with
82 % within a factor of 10 and ρ = 0.75 (p=0.058).

An undisturbed sample takes Ks only from undisturbed neighbours. That follows
soil physics — intact soil keeps its macropores — but changes little in
practice, since 10,951 of the 11,017 reference layers with a measured Ks are
undisturbed: the paired error is much the same either way (0.87 against 0.86
dex for a new source, p=0.53; 0.56 against 0.56 with the source in the
reference, p=0.59). Telling the classifier the type barely moves texture for
undisturbed soils (+0.2 pp, p=1.00, new source). The reference
holds only 45 disturbed layers with a measured Ks (ETH literature
compilations, UNSODA and AfSPDB), so a **disturbed** sample takes Ks from all
soils, and the tool says so.

### Local covariates (depth, bulk density)

A user usually knows the sensor depth and often the bulk density of a core
from the same spot. Both were tested as extra inputs to the classifier, and
as extra matching dimensions for the neighbours (`verify_covariates.py`;
1,181 targets from 36 laboratories, paired against the default):

| added input | source already in the reference | new source |
|---|---|---|
| depth | +2.5 pp (p=0.039) | +0.7 pp (p=0.55) |
| bulk density | +2.6 pp (p=0.010) | +0.0 pp (p=1.00) |
| depth + bulk density | **+3.6 pp (p=0.003)**; European soils −0.5 pp (p=0.91) | +0.8 pp (p=0.45) |

Exact class, as inputs to the classifier. With the source in the reference
the texture group moves the same way (+2.4 pp with both, p=0.051), and at the
standard leave-one-soil-out level the pair adds +3.8 pp (p=0.003). **They help
when the user's data source is in the reference; for a new source the effect
is small and not stable**: +0.8 pp here, −0.1, +0.9, +2.5 (p=0.032) and
−1.0 pp (p=0.39) on earlier draws of test soils, and the texture group does not move
(+0.2 pp, p=0.94). The gain with the source in the reference ranges from
+2.8 to +4.5 pp across draws, and for European soils it has come and gone
(−0.5 pp here, +3.8 pp on an earlier draw).
They are optional inputs (`--depth`, `--bulk-density`), worth supplying,
most of all once your own verified data have been added.

Part of the gain is a source signature: depth places a sample in its profile
and bulk density describes its packing, but both also carry how a source
samples and prepares its soils, so the classifier leans towards the classes
that source contains most. Testing the Zanjanrood soils with their own data
in the reference, depth and bulk density lifted clay loam (the commonest
class there) from 56 % to 83 % correct but cut loam from 25 % to 4 %
(`verify_external.py`). Expect better accuracy overall and worse on the
classes your verified data rarely contain.

As an extra matching dimension for the neighbours, depth changes neither
texture nor Ks. Bulk density there lowers the Ks error slightly with the
source in the reference (0.41 → 0.39 dex, p=0.07; leave-one-soil-out
0.40 → 0.37 dex, p=0.01) and for a new source (1.10 → 0.97 dex, p=0.003).
The tool therefore matches the neighbours on bulk density whenever you
supply it; depth reaches the classifier only.

### Volcanic and andic soils

Volcanic soils are harder than others: from a new source, 21 % exact and
42 % texture group against 28 % and 56 % for other soils
(`verify_volcanic.py`, 500 volcanic soils and 1,157 others, three random
draws). Andic soils disperse poorly and hold much water for their clay, so
their curves read as coarser or siltier than their measured texture.

Stating the soil's andic properties (`--andic`) is tested through the whole
pipeline (`verify_andic.py`; each target's source held out, the rest of
the reference kept):

| | exact class | texture group | fractions MAE | Ks |
|---|---|---|---|---|
| 341 andic soils, new source | +14.1 pp (p<0.001) | +20.8 pp (p<0.001) | −1.9 points (p<0.001) | unchanged |
| 480 other soils stating "no" | +2.5 pp (p=0.18) | −0.8 pp (p=0.72) | −0.1 | unchanged |
| Canary Islands andic soils, new source | −10.6 pp (p=0.039) | +21.2 pp (p=0.009) | +0.1 (p=0.69) | unchanged |

The size of the gain depends on the rest of the reference: before the Yellow
River soils joined it was +8.8 pp exact and +12.0 pp group for the same 341
andic soils. With the source in the reference it is +18.2 pp exact (+7 to
+11 pp on earlier references) and −2.4 points fractions. For other soils,
stating "no" with the source in the reference has moved exact class by −6.2
to +2.1 pp across draws (−2.9 pp, p=0.16, here): no consistent effect, so
leaving the flag out costs nothing. The classifier takes the flag as a feature; the
neighbours add one standard unit of distance between soils that differ in
it, which moves the class vote and the fractions but not Ks — when this was
tested the only andic soils with a measured Ks came from one source
(Campania), so matching Ks on it helped that source and misled every other
(`verify_andic_knn.py`); the Canary Islands soils now add a second, not yet
retested. The gain also depends on andic soils from more than
one source: with KSSL and Campania held out together (278 of the 341), the
flag cost 11 points. The Canary soils' exact class is scored against a
hexametaphosphate texture that under-disperses allophane.

Each was tested with paired designs; the scripts are in `dev/`.

| change | result | why not adopted |
|---|---|---|
| Merge UNSODA 2.0 or sDB | +0.7 pp (p=0.052); +0.3 pp (p=0.83) | mostly already inside GSHP |
| Organic carbon as a feature | +1.3 pp (p=0.33) | measured for only 42 % of the reference |
| Match on curve shape instead of parameters | −0.6 pp; −1.3 pp whitened | the eight curve points carry ~1.6 independent dimensions |
| Free m from n (five parameters) | −3.0 pp kNN (p=0.007) | a fifth axis thins the neighbourhood |
| Temper the class prior | +0.7 pp (p=0.33) | not significant; the uniform prior is kept |
| Weight neighbours by fit quality | +0.0 pp | GSHP fits are mostly well constrained |
| More European data (EU-HYDI) for an unseen European laboratory | +2.6 pp (p=0.19) | kept for coverage; the gap is laboratory protocol, not geography |
| Volcanic parent material as a classifier feature | +7 to +11 pp with the source in the reference; texture group −7 to −35 pp for non-andic volcanic soils from a new source | acts as a source label; the andic flag carries the useful part |
| Andic properties in the Ks neighbour search | Ks +0.1 to +1.0 dex worse for a new source | one source held every andic Ks (before the Canary Islands joined) |
| Stronger andic matching (2 units, or andic-only neighbours) | fractions −2 to −3 points for andic soils, but Canary exact −14 pp (p=0.012) | the one-unit penalty keeps most of the gain without it |

What *has* helped: the gradient-boosted classifier out-of-laboratory, the
Hohenbrink cores (Ks band coverage on those cores from 51 % to 65 %),
KSSL's class mix (1.8 times GSHP's silty clay loam; its accuracy gain alone
was not significant), and a laboratory's own data — having the
target's laboratory in the reference is worth +5 points, which is what
`add_local_data.py` is for.

## Files

**The tool**

- `swcc_texture.py` — fit, inference and command line
- `add_local_data.py` — add your own lab-verified curves to the reference
- `data/*_reference.csv` — the reference tables (EU-HYDI's built locally only)
- `Testing/` — example curves: `<i>_<code>_CnP.csv` (Carsel & Parrish class
  means), `<i>_<code>_GHSP.csv` (real GSHP soils; they are in the reference,
  so the tool finds each one itself) and `13_SiL_UNSODA.csv` (a real
  undisturbed silt loam that is not in the default reference)

**Building the reference** — each script documents where to obtain its raw
input: `prepare_gshp.py`, `prepare_kssl.py`, `prepare_hohenbrink.py`,
`prepare_babaeian_zanjanrood.py`, `prepare_babaeian_az.py`,
`prepare_armas.py`, `prepare_tong.py`, `prepare_euhydi.py`, `prepare_willard.py`,
`prepare_unsoda.py`, `prepare_sdb.py`; `prepare_volcanic.py` (the volcanic and andic fields);
`free_m.py` (an unconstrained companion fit written alongside the Mualem
one).

**Verification** — `verify_holdout.py` (headline accuracy at each hold-out
level), `verify_hybrid.py`, `verify_ceiling.py`,
`verify_source_blocked.py` (per-laboratory leakage), `verify_sample_type.py`,
`verify_covariates.py`, `verify_external.py` (an outside dataset before it
joins the reference), `verify_volcanic.py`, `verify_andic_knn.py` and
`verify_andic.py` (volcanic and andic soils), the benchmarks `verify_carsel_parrish.py`,
`verify_rosetta.py`, `verify_gshp.py`, `verify_groups.py`, and the shared
`ksat_metrics.py` and `verify_common.py`.

**Figures** — `figures/workflow.html` (how the reference is built, how a
prediction is made, what the validation proves) with its SVGs, and
`figures/make_coverage_map.py` and `figures/make_class_distribution.py`.

**`dev/`** — superseded experiments, kept for the record (see `dev/README.md`).

## References

- Armas Espinel, S. (2013). Contribución al estudio de las propiedades
  físicas de suelos ándicos de las Islas Canarias. PhD thesis, Universidad de
  La Laguna, Serie Tesis Doctorales, Ciencias y Tecnologías 42.
- Armas-Espinel, S., Hernández-Moreno, J.M., Muñoz-Carpena, R. and Regalado,
  C.M. (2003). Physical properties of "sorriba"-cultivated volcanic soils from
  Tenerife in relation to andic diagnostic parameters. *Geoderma* 117:297–311.
- Babaeian, E., Homaee, M., Vereecken, H., Montzka, C., Norouzi, A.A. and
  van Genuchten, M.Th. (2015). A comparative study of multiple approaches for
  predicting the soil–water retention curve: hyperspectral information vs.
  basic soil properties. *Soil Science Society of America Journal*
  79:1043–1058. doi:10.2136/sssaj2014.09.0355.
- Carsel, R.F. and Parrish, R.S. (1988). Developing joint probability
  distributions of soil water retention characteristics. *Water Resources
  Research* 24(5):755–769. doi:10.1029/WR024i005p00755.
- Cover, T. and Hart, P. (1967). Nearest neighbor pattern classification.
  *IEEE Transactions on Information Theory* 13(1):21–27.
- Global Volcanism Program (2026). Volcanoes of the World (v. 5.4.0,
  7 Aug 2026). Smithsonian Institution, compiled by E. Venzke.
  doi:10.5479/si.GVP.VOTW5-2026.5.4.
- Gupta, S., Papritz, A., Lehmann, P., Hengl, T., Bonetti, S. and Or, D.
  (2022). Global Soil Hydraulic Properties dataset based on legacy site
  observations and robust parameterization. *Scientific Data* 9:444.
  doi:10.1038/s41597-022-01481-5. Data: doi:10.5281/zenodo.6640246 (CC BY 4.0).
- Hohenbrink, T.L., Jackisch, C., Durner, W., Germer, K., Iden, S.C.,
  Kreiselmeier, J., Leuther, F., Metzger, J.C., Naseri, M. and Peters, A.
  (2023). Soil water retention and hydraulic conductivity measured in a wide
  saturation range. *Earth System Science Data*.
- Nemes, A., Schaap, M.G., Leij, F.J. and Wösten, J.H.M. (2001). Description of
  the unsaturated soil hydraulic database UNSODA version 2.0. *Journal of
  Hydrology* 251:151–162.
- Leenaars, J.G.B., van Oostrum, A.J.M. and Ruiperez Gonzalez, M. (2014).
  Africa Soil Profiles Database, version 1.2. ISRIC Report 2014/01,
  ISRIC – World Soil Information, Wageningen.
- Poggio, L., de Sousa, L.M., Batjes, N.H., Heuvelink, G.B.M., Kempen, B.,
  Ribeiro, E. and Rossiter, D. (2021). SoilGrids 2.0: producing soil
  information for the globe with quantified spatial uncertainty. *SOIL*
  7:217–240. doi:10.5194/soil-7-217-2021.
- Schaap, M.G., Leij, F.J. and van Genuchten, M.Th. (2001). ROSETTA: a computer
  program for estimating soil hydraulic parameters with hierarchical
  pedotransfer functions. *Journal of Hydrology* 251:163–176.
- Soil Survey Staff, NRCS, USDA. *National Cooperative Soil Survey Soil
  Characterization Database* (Kellogg Soil Survey Laboratory).
  https://ncsslabdatamart.sc.egov.usda.gov/ (accessed 2026-08-09). Public
  domain.
- van Genuchten, M.Th. (1980). A closed-form equation for predicting the
  hydraulic conductivity of unsaturated soils. *Soil Science Society of
  America Journal* 44:892–898.
- Tong, Y., Wang, Y., Zhou, J., Guo, X., Wang, T., Xu, Y., Sun, H., Zhang, P.,
  Li, Z. and Lauerwald, R. (2024). Dataset of soil hydraulic parameters in
  the Yellow River Basin [dataset]. PANGAEA, doi:10.1594/PANGAEA.965004
  (CC BY 4.0).
- Vereecken, H., Van Looy, K., Weynants, M. and Javaux, M. (2017). Soil
  retention and conductivity curve data base sDB. PANGAEA,
  doi:10.1594/PANGAEA.879233 (CC-BY-3.0).
- Vingiani, S., Buonanno, M., Coraggio, S. et al. (2018). Soils of the Aversa
  plain (southern Italy). *Journal of Maps* 14:312–320.
  doi:10.1080/17445647.2018.1458338.
- Weynants, M. et al. (2013). European HYdropedological Data Inventory
  (EU-HYDI). EUR 26053 EN, Publications Office of the European Union.
  doi:10.2788/5936. Restricted to consortium members.
- Willard, L.L., Muñoz-Carpena, R., Venort, T., Gitonga, J., Maltais-Landry,
  G., Caylor, K.K. and Palm, C.A. A high-resolution comprehensive database to
  evaluate agricultural impacts at edge of savanna in Kenya. Under review,
  *Scientific Data*.

## License

Released under the Creative Commons Attribution 4.0 International license
(CC BY 4.0); full text in `LICENSE`. Suggested attribution: "R. Muñoz-Carpena,
py_HB_texture, https://github.com/carpena1/py_HB_texture".

The bundled GSHP-derived table (`data/gshp_reference.csv`) is itself CC BY 4.0
— cite Gupta et al. (2022) and Zenodo record 6640246 when reusing it. The sDB
table is CC-BY-3.0 (Vereecken et al. 2017). The Zanjanrood table
(`data/babaeian_zanjanrood_reference.csv`) holds E. Babaeian's own measurements,
distributed with the data owner's permission; cite Babaeian et al. (2015). EU-HYDI data and anything derived
from it at sample level are not part of this repository, nor are the
unpublished Laikipia data.
