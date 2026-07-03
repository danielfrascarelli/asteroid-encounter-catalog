# A catalogue of real 3D close encounters between asteroids during the Gaia DR3 window — draft

> **Estado:** 🟡 DRAFT — 2026-07-03
> Draft del dataset paper (frente P4 de [`planning/PUBLISH_PUSH_PLAN.md`]).
> Abstract, §1–§3 y §5–§7 redactadas; §4 llena con la minería P1
> (`docs/notable_encounters.md`). Marcadores «[TODO: …]» señalan cifras
> pendientes de corridas en curso (Metis F4/F3).

## Título tentativo

A systematic catalogue of real three-dimensional close encounters between
numbered asteroids during the Gaia DR3 observation window, with a measured
completeness budget.

---

## Abstract

We present a systematic catalogue of real three-dimensional close encounters
between numbered asteroids during the Gaia DR3 observation window
(2014 July 25 – 2017 May 28). Orbits are propagated from a single frozen MPCORB
snapshot whose osculating epoch sits at the centre of the window, and encounters
are detected by a per-timestep KD-tree spatial query, refined on a dense temporal
sub-grid around each apparent minimum. The catalogue records the minimum physical
3D separation of each pair — a genuine spatial approach, not an apparent
sky-plane co-location — for 72,236,904 pairs that came within 0.05 AU, together
with the encounter epoch, relative velocity, observing geometry, and estimated
physical properties. A hybrid variant carries full N-body (rebound/ASSIST)
refinement for the 8.7 million pairs in the dynamically fragile subset
(low perihelion or high eccentricity). The distinguishing feature of this work
is a *measured* completeness budget rather than an assumed one: the Kepler-to-N-body
refinement error (median 12 μAU, p99 2.5 mAU), the false-negative rate induced by
the 0.05 AU threshold (0.70 %, 95 % CI 0.59–0.83 %, symmetric and
scatter-dominated), and the recall of the orbital prefilter on the adverse
high-e/high-i tail (76 %) are each quantified from dedicated re-refinement
experiments. We show the catalogue is a direct target-selection input for
asteroid mass determination: applying a joint orbit+mass least-squares engine to
Gaia FPR astrometry recovers all four calibrator masses (Ceres, Vesta, Pallas,
Hygiea) at |z| < 3 and yields sixteen perturber determinations consistent with
recent independent work.

---

## 1. Introduction

Close encounters between asteroids are of interest on several fronts. A slow,
close approach between a massive perturber and a small test body deflects the
latter's orbit measurably, and such events are the classical route to asteroid
mass determination (Michalak 2000; Goffin 2014; Fuentes-Muñoz et al. 2024, 2025).
Encounters between members of the same dynamical family probe the collisional
history of the belt, and the statistical distribution of mutual approach
distances and velocities is itself a dynamical observable.

Existing work on encounter-based mass determination — Goffin (2014),
Fuentes-Muñoz et al. (2024, 2025) — starts from *hand-selected* individual
encounters chosen for their favourable geometry, typically a single large
perturber against a short list of well-observed test bodies. There is, to our
knowledge, **no systematic catalogue of real 3D asteroid–asteroid encounters
whose own completeness has been characterised**. That is the gap this work fills:
we detect encounters exhaustively over the full numbered population in the Gaia
DR3 window, and — the point that separates this catalogue from a bulk
propagation exercise — we *measure* what the detection pipeline misses rather
than asserting completeness.

A methodological distinction is central. Two asteroids can appear arbitrarily
close in right ascension and declination as seen from Earth while remaining
astronomical units apart along the line of sight. This catalogue detects the
former case only when it is also the latter: it records the **minimum physical
separation in three dimensions** between the two bodies' heliocentric positions,
independent of any observer. It is therefore a catalogue of genuine spatial
approaches, not of sky-plane co-localisations. Light-time correction — essential
when comparing model positions to Gaia's astrometry — is deliberately *not*
applied here, because the encounter is a purely geometric relation between two
propagated trajectories at a common dynamical time.

The remainder of this paper describes the data and detection method (§2), the
measured completeness budget that is the catalogue's distinguishing feature (§3),
the notable individual events surfaced by mining the catalogue (§4), the
application to mass determination (§5), data availability (§6), and conclusions
(§7).

---

## 2. Data and method

### 2.1 Sources and conventions

Orbital elements are taken from a single frozen MPCORB snapshot,
`MPCORB_20160217.DAT` (SHA-256 prefix `3e44e7d3…`), selected because its
osculating epoch lies at the centre of the Gaia DR3 window, minimising the
propagation distance — and hence the accumulated two-body error — for Kepler
refinement. The snapshot hash and download date are recorded in the catalogue's
provenance sidecar so that results are traceable to an exact input; MPCORB is a
living file and results are not bit-reproducible across snapshot versions. The
run is restricted to **numbered asteroids** (higher orbital quality); provisional
designations are out of scope for this release.

The temporal domain is the Gaia DR3 Solar System survey window,
2014 July 25 – 2017 May 28. All internal computation is carried out in Julian
Date on the Barycentric Dynamical Time (TDB) scale; conversions to/from the
Barycentric Coordinate Time (TCB) that Gaia reports, and to UTC, are performed
only at input/output interfaces. Distances are in AU, angles in radians, and
velocities in AU/day internally, converted to user-facing units (km, degrees,
km/s) only for output.

### 2.2 Detection

The naive cost is prohibitive: ~150,000 bodies over a multi-year window at
hour-level cadence implies of order 10¹⁰ pairwise distances per step. Detection
is therefore layered.

An **orbital prefilter** is available to prune pairs that cannot physically
approach — the shipped heuristic keeps a pair only if |Δa| ≤ 0.5 AU and
|Δi| ≤ 30°. This filter is a deterministic mask on osculating elements. We note
here, and quantify in §3.3, that at main-belt scale (≫ 5000 bodies) the pair-list
prefilter is not built; the frozen run therefore applied only the KD-tree spatial
query, and the |Δa| cut did *not* shape this catalogue (provenance
`prefilter.effective = "skipped_large_n"`).

The core of the detector is a **per-timestep 3D KD-tree**. At each coarse step
the heliocentric ecliptic positions of all bodies are inserted into a KD-tree
(O(N log N) construction), and a radius query at the 0.05 AU threshold returns
neighbour candidates (O(log N) per query), avoiding any O(N²) distance matrix.
The coarse cadence is Δt = 12 h.

Each coarse candidate is then **refined on a dense temporal sub-grid**: a
±2 h window around the apparent minimum is sampled at 120 s, and the minimum
separation, encounter epoch, and relative velocity are extracted from the
refined track. Time blocks are independent, so the scan parallelises across
processes.

**Propagation** is two-tiered. The coarse KD-tree scan is driven by an N-body
rebound integration (WHFast; Sun + Jupiter + Saturn; dt = 1 h). The sub-grid
refinement that produces the *reported* minimum distance runs through an
analytic two-body Kepler propagator. Over the short Gaia DR3 arc (~3 yr from the
central epoch) two-body propagation is adequate for the bulk of the belt, and its
residual error is measured, not assumed (§3.1). For the dynamically fragile
subset — pairs with perihelion `q_min < 1.8 AU` or eccentricity `e_max > 0.3`,
8,728,509 pairs (12.08 % of the catalogue) — a **hybrid** catalogue replaces the
Kepler distances with full N-body refinement (rebound IAS15; Sun + 8 planets + 4
major asteroids; ±12 h around the Kepler minimum; 60 s sampling). The hybrid
catalogue carries `refinement_method ∈ {kepler, nbody}` per row together with the
Kepler and N-body values and their residuals, so users needing sub-mAU geometry
on high-e / NEA-crossing pairs can cite it directly. Universal N-body refinement
of all 72 M pairs (~400 CPU-days) is explicitly out of scope: the median error
outside the fragile subset is already below 1 mAU (§3.1), so the marginal gain
does not justify the cost.

### 2.3 Characterisation

For each detected pair the catalogue records the minimum separation and its exact
epoch (from the refined sub-grid), the relative velocity at closest approach, and
the geometry of the event relative to Gaia — solar elongation and an
observability flag — computed on the characterised catalogue. Of the 72,236,904
encounters, 13,640,870 (18.9 %) fall in a Gaia-observable configuration. Physical
properties (absolute magnitude H, and diameter and taxonomic class estimated from
H with a class-dependent albedo when no direct measurement exists) are attached
to both bodies. Diameters derived this way are order-of-magnitude estimates, not
measurements, and any size-based selection inherits that uncertainty.

---

## 3. Completeness budget (the distinguishing feature)

The value of this catalogue is not that it detects encounters — that is a bulk
propagation — but that it states, with measured numbers, what it *misses* and
how large the Kepler-vs-N-body distortion is. Three quantities are measured;
each is reproducible from a dedicated script and documented in a companion note.

### 3.1 Kepler vs N-body refinement error

The reported distances in the frozen catalogue come from the two-body refiner.
We measured the resulting error against full N-body refinement in two stages
([`docs/kepler_refine_error_report.md`]):

- **Stage A** (964 stratified pairs re-refined under full N-body): median
  |Δdist| = **12 μAU**, p95 = **678 μAU**, p99 = **2.5 mAU**, max = **11.3 mAU**.
- **Stage B** (production N-body refinement of the entire 8.73 M fragile subset):
  p99 |Δdist| = **1.99 mAU**, max = **15.2 mAU**; 0 failed integrations, 0
  unconverged, maximum energy drift 5.6 × 10⁻¹⁴.

For scale, 1 mAU ≈ 1.5 × 10⁵ km, and the detection threshold is 0.05 AU. The
median distortion is thus ~4 orders of magnitude below the threshold, and even
the Stage B worst case is ~30 % of a threshold — confined to the fragile subset,
which the hybrid catalogue refines. A separate check of the truncated perturber
set (the scan used only Sun + Jupiter + Saturn) shows that adding the remaining
planets shifts closest-approach distances by a median 1.3 μAU (p99 67 μAU, max
80 μAU), i.e. ~100× below the Kepler two-body error — the perturber truncation is
not the dominant error term ([`docs/nbody_perturber_ceiling.md`]).

### 3.2 Threshold censoring (measured false-negative rate)

Because the catalogue contains only pairs whose *Kepler* distance is below
0.05 AU, a subtle bias arises at the threshold. When the fragile subset is
re-refined to N-body, 25,283 pairs (0.29 %) cross the threshold *upward*
(Kepler < 0.05, N-body ≥ 0.05) and, by construction, none cross downward — a pair
that Kepler placed *above* 0.05 AU was never written, so its downward crossing
**cannot be observed**. The apparent one-sided "no false negatives" is therefore
censoring, not measurement ([`docs/kepler_threshold_bias_paper.md`]).

We measured the missing side directly. Re-running Kepler detection at a widened
0.06 AU threshold over 10,000 numbered bodies and isolating the 17,469 pairs with
Kepler distance in [0.05, 0.06) AU — exactly the band a 0.05 AU catalogue
discards — then re-refining each under full N-body, **0.70 %** [95 % CI
0.59–0.83 %] cross *downward* below 0.05 AU (a conservative floor of ~0.42 %
after excluding near-boundary cases). The residual Δdist (N − K) has median
−3 × 10⁻⁷ AU and σ ≈ 3.7 × 10⁻⁴ AU — a symmetric, scatter-dominated distribution,
not a directional bias. In other words the crossing is a threshold selection
effect (Eddington/Malmquist family): scatter of amplitude ~1 σ ≈ 4 × 10⁻⁴ AU
pushes pairs just inside the threshold to just outside, and their mirror images
(just outside to just inside) fall outside the catalogue. The upward rate
(~1.5 % in the boundary bin [0.045, 0.050)) and the downward rate (~0.4–0.7 % in
[0.05, 0.06)) are consistent with a single symmetric process. Crossings
concentrate in low-q / high-e / fast orbits, exactly where two-body propagation
is weakest — which independently validates the `q_min < 1.8 ∨ e_max > 0.3`
subset criterion.

Extrapolated catalog-wide (order-of-magnitude, the [0.05, 0.06) band scaling as
N²), this implies of order **10⁵ real < 0.05 AU encounters censored** by the
Kepler cut (~1.5–2.5 × 10⁵), separate from and additional to the prefilter
recall deficit of §3.3. This is why the word "complete" is never used for this
catalogue.

### 3.3 Prefilter recall on the adverse tail

The shipped |Δa| ≤ 0.5 AU prefilter is blind to eccentricity: a high-e orbit
reaches well inside the belt at perihelion and well outside at aphelion, so it
can physically cross an orbit whose semimajor axis differs by far more than
0.5 AU. We measured the damage this would cause on the adverse subset (numbered,
a ∈ [1.5, 4.0], e > 0.3 ∨ i > 15°; 52,411 bodies), exploiting the fact that the
prefilter is a pure deterministic mask: running detection without the prefilter
and intersecting with the mask is byte-identical to running with it
([`docs/prefilter_recall.md`]). Ground truth is 606,393 real adverse–adverse
encounters < 0.05 AU; the prefilter keeps 463,164, for a recall of **76.38 %**
[95 % CI 76.27–76.49 %]. **143,229** encounters would be missing from a
prefiltered catalogue (a lower bound — adverse–normal pairs were not measured),
and 99.6 % of that loss is the |Δa| cut alone. Recall degrades smoothly with
eccentricity, from ~99.9 % below e = 0.1 to ~39 % above e = 0.7.

Two points bound the impact. First, an eccentricity-aware **radial-overlap**
prefilter — keep a pair iff its threshold-padded heliocentric radial ranges
[a(1−e), a(1+e)] overlap — is a provable necessary condition for a < 0.05 AU
encounter, is just as cheap, and recovers **100 %** recall on the adverse subset;
it is the recommended replacement for any future run. Second, and specific to
*this* freeze: because the run was at main-belt scale the pair-list prefilter was
skipped entirely, so the |Δa| cut did not actually drop these pairs from the
frozen catalogue. The 76 % figure measures the damage the prefilter *would* cause
if applied at small N, not damage suffered here. Completeness of this freeze is
instead bounded by the 12 h KD-tree coarse cadence (fast minima can fall between
grid samples) and by the threshold censoring of §3.2. For the dynamically cold
bulk of the belt (low e, small Δa) recall is ~99.9 %; the incompleteness is a
property of the high-e/high-i tail, which any science touching NEA-crossing or
high-e pairs must cite.

---

## 4. Notable events

Mining the characterised catalogue (`encounters_characterized_full.parquet`;
`scripts/bench/mine_notable_encounters.py`, [`docs/notable_encounters.md`])
surfaces individual events across four axes. We stress at the outset that the
value of these tables is statistical and archival rather than a single headline
discovery: there is **no new large–large (D ≳ 100 km) encounter** — the two such
events in the catalogue are both already known — and every candidate below is
under the frozen two-body assumptions and requires N-body revalidation before
publication. Distances are physical 3D separations at closest approach; diameters
are H-derived estimates.

### 4.1 Large–large encounters

Only two encounters have both bodies with estimated D ≳ 100 km, and both involve
already-studied bodies:

| body 1 | body 2 | date (UTC) | dist (km) | v_rel (km/s) | D₁ (km) | D₂ (km) |
|---|---|---|---:|---:|---:|---:|
| (1) Ceres | (57) Mnemosyne | 2017-01-10 | 6.5 × 10⁶ | 7.25 | 763 | 139 |
| (7) Iris | (44) Nysa | 2014-08-13 | 7.3 × 10⁶ | 5.35 | 281 | 139 |

The population of both-D ≳ 50 km encounters is richer (top by minimum distance):

| body 1 | body 2 | date (UTC) | dist (km) | v_rel (km/s) | D₁ (km) | D₂ (km) |
|---|---|---|---:|---:|---:|---:|
| (305) Gordonia | (830) Petropolitana | 2014-09-18 | 9.88 × 10⁵ | 3.28 | 62.6 | 53.8 |
| (426) Hippo | (788) Hohensteina | 2015-04-17 | 1.15 × 10⁶ | 8.82 | 73.5 | 64.6 |
| (145) Adeona | (780) Armenia | 2016-11-09 | 1.50 × 10⁶ | 6.34 | 84.0 | 56.3 |
| (739) Mandeville | (415) Palatia | 2017-04-05 | 1.70 × 10⁶ | 8.36 | 70.9 | 51.1 |
| (675) Ludmilla | (80) Sappho | 2016-02-15 | 1.71 × 10⁶ | 2.97 | 93.0 | 90.0 |
| (51) Nemausa | (91) Aegina | 2016-07-13 | 2.08 × 10⁶ | 4.42 | 120 | 60.6 |
| (46) Hestia | (482) Petrina | 2016-07-12 | 2.46 × 10⁶ | 4.23 | 75.6 | 60.6 |
| (758) Mancunia | (740) Cantabia | 2016-08-13 | 2.63 × 10⁶ | 3.70 | 82.9 | 53.8 |
| (739) Mandeville | (754) Malabar | 2016-03-14 | 2.72 × 10⁶ | 5.26 | 70.9 | 51.6 |
| (200) Dynamene | (1304) Arosa | 2014-12-25 | 2.85 × 10⁶ | 7.88 | 79.2 | 51.3 |

The closest both-D ≳ 50 km event, (305) Gordonia × (830) Petropolitana at
9.9 × 10⁵ km, involves neither of the sixteen studied perturbers and is a natural
target for N-body revalidation.

### 4.2 Extreme events

The absolute closest approaches in the catalogue are between small bodies, at
separations of ~1000–7000 km:

| body 1 | body 2 | date (UTC) | dist (km) | v_rel (km/s) | D₁ (km) | D₂ (km) |
|---|---|---|---:|---:|---:|---:|
| (153222) 2000 YD43 | (238587) 2004 YX3 | 2016-03-11 | 1094 | 6.90 | 2.82 | 1.78 |
| (15072) Landolt | (387599) 2001 XF180 | 2014-12-08 | 1760 | 5.87 | 3.24 | 1.62 |
| (270730) 2002 QE130 | (366918) 2005 UC211 | 2016-05-29 | 2291 | 3.10 | 1.62 | 1.18 |
| (117065) 2004 KD9 | (439086) 2011 QP5 | 2015-07-27 | 2535 | 8.25 | 3.72 | 1.86 |
| (52249) 1981 EK21 | (408138) 2013 CL75 | 2016-11-17 | 2590 | 3.16 | 3.39 | 1.41 |

At these separations two-body propagation is most fragile, so these are precisely
the events most in need of N-body revalidation; those that survive are the
geometrically most striking events in the dataset.

The slowest encounters (relative velocity down to ~6 m/s, with at least one body
D ≳ 5 km) are the natural targets for future mass determination, since a low
v_rel prolongs the gravitational interaction and enlarges the deflection:

| body 1 | body 2 | date (UTC) | dist (km) | v_rel (km/s) | D₁ (km) | D₂ (km) |
|---|---|---|---:|---:|---:|---:|
| (3791) Marci | (110964) 2001 UW170 | 2017-01-18 | 5.85 × 10⁶ | 0.006 | 13.5 | 2.46 |
| (47042) 1998 WP3 | (99861) Tscharnuter | 2014-07-27 | 6.03 × 10⁶ | 0.015 | 5.63 | 3.24 |
| (7203) Sigeki | (73664) 1981 EE34 | 2015-02-28 | 7.07 × 10⁶ | 0.031 | 8.14 | 1.48 |
| (3384) Daliya | (309262) 2007 RZ93 | 2016-12-20 | 4.44 × 10⁶ | 0.031 | 6.17 | 1.23 |
| (7479) 1994 EC1 | (14954) 1996 DL | 2017-05-12 | 7.10 × 10⁶ | 0.032 | 6.46 | 3.09 |

The combined "large + slow + close" set (one body D ≳ 50 km, v_rel ≤ 1 km/s,
dist ≤ 10⁶ km) is of the highest physical interest, being the configuration that
maximises measurable deflection:

| body 1 | body 2 | date (UTC) | dist (km) | v_rel (km/s) | D₁ (km) | D₂ (km) |
|---|---|---|---:|---:|---:|---:|
| (135) Hertha | (281202) 2007 FE47 | 2016-09-18 | 4.56 × 10⁵ | 0.579 | 80.3 | 1.02 |
| (371) Bohemia | (14711) 2000 CG36 | 2016-05-14 | 5.92 × 10⁵ | 0.982 | 64.0 | 5.90 |
| (83) Beatrix | (170676) 2003 YE180 | 2014-09-18 | 6.10 × 10⁵ | 0.684 | 65.8 | 2.14 |
| (678) Fredegundis | (200373) 2000 QH89 | 2015-03-13 | 6.27 × 10⁵ | 0.964 | 55.8 | 2.04 |
| (9) Metis | (345254) 2005 UT477 | 2016-10-30 | 7.14 × 10⁵ | 0.979 | 197 | 1.23 |
| (639) Latona | (10793) Quito | 2014-12-18 | 7.27 × 10⁵ | 0.604 | 81.4 | 10.7 |
| (16) Psyche | (139485) 2001 PT14 | 2015-01-19 | 9.66 × 10⁵ | 0.869 | 235 | 2.57 |
| (7) Iris | (153539) 2001 SD101 | 2015-08-17 | 9.97 × 10⁵ | 0.806 | 281 | 2.35 |

### 4.3 Same-region pairs (family proxy)

Pairs whose osculating elements are mutually close (Δa/a ≤ 1 %, Δe ≤ 0.02,
Δi ≤ 1°) *and* that had a physical close approach are candidate co-orbital /
same-family associations. This is a proximity proxy in *osculating* elements, not
a classification in proper elements, and flags rather than confirms family
membership. The closest such pair, (200764) 2001 XP3 × (131597) 2001 XT2, passed
within 49,511 km on 2014-08-13 with nearly identical (a, e, i) —
a ≈ 2.29 AU, e ≈ 0.19, i ≈ 25° for both — a plausible genetic pair worth
follow-up in proper-element space.

### 4.4 Candidate perturbers beyond the studied sixteen

Ranking large bodies (D ≳ 100 km) *not* among the sixteen studied perturbers by
their number of *useful* encounters (v_rel ≤ 3 km/s and dist ≤ 3 × 10⁶ km) yields
the best prospects for a mass determination the frozen mass runs have not yet
attempted:

| body | D (km) | class | useful events | best dist (km) | min v_rel (km/s) |
|---|---:|---|---:|---:|---:|
| (9) Metis | 197 | MBA | 37 | 7.14 × 10⁵ | 0.378 |
| (30) Urania | 109 | — | 36 | 4.40 × 10⁵ | 0.817 |
| (40) Harmonia | 141 | — | 36 | 6.58 × 10⁵ | 1.028 |
| (19) Fortuna | 133 | MBA | 32 | 7.43 × 10⁵ | 0.678 |
| (21) Lutetia | 120 | — | 30 | 8.18 × 10⁵ | 0.609 |
| (64) Angelina | 104 | MBA | 30 | 8.73 × 10⁵ | 0.992 |
| (44) Nysa | 140 | — | 30 | 8.95 × 10⁵ | 0.961 |
| (20) Massalia | 178 | — | 28 | 8.78 × 10⁵ | 0.514 |
| (29) Amphitrite | 240 | — | 27 | 3.17 × 10⁵ | 0.643 |
| (27) Euterpe | 141 | — | 26 | 4.63 × 10⁵ | 0.847 |

(9) Metis heads the list with 37 useful events and a best pair at v_rel = 0.38 km/s;
(29) Amphitrite is the largest (D ≈ 240 km) with a close 3.2 × 10⁵ km approach.
These candidates directly informed the F4 mass determination of §5.

### 4.5 Separating new from known

A reference list of encounters touching one of the sixteen studied perturbers
(1, 2, 3, 4, 7, 10, 15, 16, 31, 52, 65, 87, 88, 107, 511, 704) is provided in
[`docs/notable_encounters.md`] so that already-worked events can be separated from
potential discoveries. The honest verdict is that no new large–large encounter of
publication interest emerged: the belt's largest bodies are well studied, and
Fuentes-Muñoz et al. (2025) already published 231 masses covering essentially all
classical large MBAs. The novelty of this work is the dataset itself — a catalogue
with a measured completeness budget — and the methodological framework of §5, not
a single spectacular event.

---

## 5. Application to mass determination

The catalogue is the target-selection front-end of an asteroid mass-determination
engine. Perturber mass is obtained by a joint least-squares fit of the perturber's
mass and the orbits of its test bodies — those with a < 0.05 AU encounter drawn
directly from the catalogue via `--from-catalog` — over the full arc of Gaia FPR
astrometry. The force model is ASSIST (JPL DE440 ephemeris; EIH relativistic
corrections; sixteen massive asteroid perturbers); the implementation is our own
(`src/orbdet/`), with no third-party orbit-determination dependency. Full results
are in [`docs/mass_determination_results.md`] and [`docs/mass_crosscheck_jack.md`].

**Calibrators.** The four calibrator masses are recovered at |z| < 3 against the
literature. With N ≥ 20 test bodies, the fit/literature ratios of Ceres, Vesta
(DAWN) and Hygiea (Vernazza et al. 2020) fall in [0.943, 0.990] (mean bias −4 %);
Pallas, limited to 6 encounters in the catalogue, gives ratio 1.240 (z = +2.67).

| body | N | fitted mass (kg) | σ_total | ratio fit/lit | z |
|---|---:|---:|---:|---:|---:|
| Ceres | 28 | 8.96 × 10²⁰ | 4.6 % | 0.955 | −1.01 |
| Vesta | 28 | 2.44 × 10²⁰ | 4.6 % | 0.943 | −1.30 |
| Hygiea | 20 | 8.22 × 10¹⁹ | 5.7 % | 0.990 | −0.13 |
| Pallas | 6 | 2.54 × 10²⁰ | 7.0 % | 1.240 | +2.67 |

**FOV-block covariance.** Each Gaia focal-plane crossing yields ~7 CCD
observations (~5 s apart) with correlated residuals. Treating them as independent
underestimates σ(mass) by a factor 1.66 (effective N ≈ 0.36 N). The measured
intra-crossing correlation is ICC = 0.32; we whiten each crossing with a
block covariance `C = diag(σ_AL²) + s_c²·11ᵀ`, the correlated floor `s_c`
calibrated by bisection to χ²_red = 1.

**Sixteen perturbers.** The engine determines all sixteen. Only (16) Psyche
reaches calibrator-grade precision outside the calibrators: M = 2.43 × 10¹⁹ kg,
σ_stat = 3.3 %, N = 36, χ²_red = 0.99, ratio 1.020 vs DE441 and 1.014
(z = +0.25) vs the independent Fuentes-Muñoz et al. (2025) FPR fit. Six perturbers
whose deflection sits below the per-encounter astrometric noise return ratios in
[0.39, 0.72]: this is the mass↔orbit degeneracy (the fit reproduces the
astrometry with a smaller mass and a re-adjusted orbit — regression toward zero
mass), and for these the formal Fisher σ underestimates the error.

**External σ by jackknife.** Because the formal σ scales as 1/√N and ignores the
mass↔orbit regression seen only by leaving encounters out, we report a jackknife σ
that reincorporates it ([`docs/mass_crosscheck_jack.md`]). It is systematically
2–3× the formal σ. Its diagnostic value is visible in the calibrators themselves:
Pallas and Vesta — whose true masses are known — sit at |z| > 3 under the formal σ
but within |z| < 3 under the jackknife σ, showing the formal bar is genuinely too
tight rather than the jackknife being permissive. Cross-checking the ten
*measured* (identifiable, snr_jack ≥ 3) perturbers against Fuentes-Muñoz (2025),
**10/10 fall within |z| < 3 under the jackknife σ, versus 5/10 under the formal σ**
(6/6 vs 3/6 for the non-calibrators, the genuinely independent comparison). The
remaining perturbers are flagged `not_identifiable` — their deflection does not
clear the per-encounter noise and their masses are reported as explicit bounds,
not measurements (a wide jackknife σ makes almost any value "consistent", which is
not the same as a detection).

**Beyond the sixteen (F4).** The engine is not restricted to the sixteen
ephemeris perturbers: each perturber is integrated as a massive rebound particle
from its own orbital elements, so a seventeenth body outside `sb441-n16` is added
as a seventeenth massive particle with no double-counting (the ASSIST `ASTEROIDS`
force is deliberately excluded). This was demonstrated on (19) Fortuna — orbit
fixed from JPL Horizons, mass free — which fits with χ²_red = 0.977 and
M = 1.13 × 10¹⁹ ± 2.2 × 10¹⁸ kg. The methodological gate is green: the engine
generalises to bodies beyond the ephemeris set. Fortuna is not itself a new mass
(it appears in Goffin 2014 and Fuentes-Muñoz 2025), and its 20 % σ makes it a
coarse determination consistent with the literature within the jackknife bar
(z = +1.25). The catalogue's top ranked out-of-sixteen candidate, (9) Metis
(37 useful encounters, §4.4), was run as a second out-of-sixteen case and fits
with χ²_red = 0.981 and M = 4.74 × 10¹⁸ ± 1.53 × 10¹⁸ kg (55 targets, jackknife
σ = 32 %), a ratio of 0.73 against Fuentes-Muñoz (2025) (6.48 × 10¹⁸ kg),
consistent within the jackknife bar (z = −1.14). Fortuna and Metis thus scatter to
either side of the literature (+30 % and −27 % in central value) but both remain
consistent once the honest external σ is used — confirming the method generalises
while underscoring that a single out-of-sixteen perturber with ~50 Gaia-FPR
targets yields only a coarse (20–30 %) determination, not a competitive mass.
Neither is a new mass: both are already in Fuentes-Muñoz (2025).

We tested whether the mean −4 % calibrator bias is an artefact of an incomplete
perturber background by enlarging it from the sixteen ephemeris bodies to
thirty-five — adding the twenty most massive asteroids of Fuentes-Muñoz (2025)
(≈1.6×10²⁰ kg spread across the belt), with masses from their Table 5 and orbits
from JPL Horizons, using the same custom-perturber machinery. The calibrator masses
shift by < 0.25 % (Ceres −0.02 %, Vesta −0.23 %, Hygiea +0.18 %), three orders of
magnitude below the deficit being probed; the systematic floor f_sys does not drop
(4.16 % → 4.26 %). Background incompleteness is therefore *not* the origin of the
−4 % bias — a far-field, dispersed perturber population averages to nearly zero
deflection over the target ensemble. The residual bias is more plausibly attributable
to per-target orbit imperfection over the short Gaia arc or unmodelled per-encounter
astrometric systematics.

---

## 6. Data availability

The primary product is the Kepler-refined candidate catalogue,
`encounters_catalog_rebound_005au.parquet` (72,236,904 rows,
SHA-256 `b0272be7…`), a typed Parquet file with a machine-readable provenance
sidecar recording the MPCORB snapshot hash, the frozen configuration, and the
pipeline parameters. The derived **hybrid** catalogue,
`encounters_catalog_hybrid_stageb.parquet`, carries per-row
`refinement_method ∈ {kepler, nbody}` and both the Kepler and N-body distances,
epochs, velocities, and their residuals; it should be cited for any science
depending on sub-mAU geometry. The characterised catalogue,
`encounters_characterized_full.parquet`, adds observability, solar elongation,
magnitudes, diameters, and taxonomic classes for the full 72 M rows. All
completeness-budget experiments are reproducible from the scripts named in §3.
Anything derived from the catalogue must cite both the freeze
(catalog SHA `b0272be7…`) and the code commit at its own generation time.

---

## 7. Conclusions

We have built and characterised a systematic catalogue of 72,236,904 real 3D
close encounters (< 0.05 AU) between numbered asteroids during the Gaia DR3
observation window. Unlike prior encounter-based work, which selects individual
events by hand, the catalogue is exhaustive over the numbered population and — its
distinguishing feature — comes with a *measured* completeness budget rather than
an assumed one: the two-body refinement error (median 12 μAU, p99 2.5 mAU), the
threshold-induced false-negative rate (0.70 %, symmetric and scatter-dominated,
implying ~10⁵ censored encounters catalog-wide), and the orbital prefilter recall
on the adverse high-e/high-i tail (76 %, with a zero-cost radial-overlap fix that
recovers 100 %). A hybrid variant provides N-body-grade geometry on the
dynamically fragile subset.

The catalogue is a working target-selection front-end for mass determination: a
joint orbit+mass engine on Gaia FPR astrometry recovers all four calibrators at
|z| < 3, determines sixteen perturbers (all ten identifiable masses consistent
with Fuentes-Muñoz 2025 at |z| < 3 under a jackknife σ), and generalises to bodies
beyond the ephemeris set. We are explicit about the limits: distances are
candidates under frozen two-body assumptions except on the N-body-refined subset,
the catalogue is not complete, and no new large–large encounter or unpublished
large-body mass emerged — the belt's massive bodies are already well surveyed. The
contribution is the dataset with its quantified budget and the methodological
framework (FOV-block covariance, identifiability criterion, external jackknife σ),
which together turn a bulk propagation into a defensibly characterised catalogue.

---

## Notas de redacción (internas)

- El gancho es la completitud MEDIDA (§3), no la lista de encuentros per se.
- Ser explícito sobre el scope: candidato bajo supuestos congelados, no completo.
- TODOs pendientes de corridas: masa de Metis (F4, §5) y f_sys extendido (F3, §5).
- Referencias a resolver en formato de journal: Michalak 2000; Goffin 2014
  (A&A 565, A56); Fuentes-Muñoz et al. 2024 (LPSC #2388), 2025 (AJ 170, 353);
  Tanga et al. 2023 (A&A 674, A12); Park et al. 2016 (DAWN Ceres); Russell et al.
  2012 (DAWN Vesta); Vernazza et al. 2020 (Hygiea).
