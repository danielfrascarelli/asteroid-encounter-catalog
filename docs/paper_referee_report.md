# Tribunal review — dataset paper on Gaia DR3 asteroid close encounters

> Referee-panel assessment of `docs/paper/aa_encounters.tex` /
> `docs/dataset_paper_draft.md`, written as a journal tribunal would judge it:
> a panel of four lenses (dynamics, statistics/completeness, catalogue/data, and an
> editor) each giving strengths, concerns, and a verdict, followed by a synthesis.
> Comparison to analogous published work is in §C and §E.
> **Date:** 2026-07-04. Judged against the merged manuscript (catalogue SHA `b0272be7…`).

---

## A. What the manuscript claims (one paragraph)

A systematic, exhaustive catalogue of **72,236,904** real 3D close-encounter pairs
(minimum heliocentric separation < 0.05 AU) among numbered asteroids over the Gaia
DR3 window (2014-07-25 – 2017-05-28), detected by a per-timestep KD-tree + temporal
sub-grid refinement, propagated two-body (Kepler) with an N-body (rebound/ASSIST)
hybrid for a dynamically fragile subset. The stated headline is not the catalogue but
its **measured completeness budget**: Kepler-vs-N-body refinement error (median
12 μAU), a threshold-induced false-negative rate (0.70 %), and a prefilter recall on
the adverse tail (76 %). A joint orbit+mass least-squares engine (ASSIST + Gaia FPR,
FOV-block covariance, jackknife σ) is applied, recovering 4 calibrator masses and 16
perturbers (10 identifiable), all consistent with recent literature.

---

## B. Referee 1 — Celestial mechanics / dynamics

**Strengths.**
- The two-body-vs-N-body error is *measured*, not asserted: median 12 μAU, p99
  2.5 mAU (Stage A), and a full production N-body re-refinement of the 8.73 M fragile
  subset (Stage B, 0 failed integrations, energy drift 5.6×10⁻¹⁴). This is the right
  way to justify a two-body candidate catalogue and it is done convincingly.
- The perturber-truncation check (Sun+Jupiter+Saturn vs full planets: median 1.3 μAU
  shift) correctly demonstrates the coarse-scan force model is not the dominant error.
- The `q_min < 1.8 ∨ e_max > 0.3` fragile-subset criterion is validated *a posteriori*
  by showing threshold crossings concentrate in low-q/high-e/fast orbits — internally
  consistent.
- Frame/time hygiene (TDB internal, TCB/UTC only at interfaces; barycentric vs
  heliocentric; light-time deliberately excluded for a geometric relation) is correct
  and explicitly stated.

**Concerns.**
- **(Raised as major → RESOLVED during this review pass.)** The submitted text said
  completeness "is bounded by the 12 h KD-tree coarse cadence (fast minima can fall
  between grid samples)" with **no number**, which for a *measured*-budget paper read
  as an unmeasured hole at the detection stage. On inspection of the pipeline this was
  a **mis-description, not a real gap**: the coarse scan *widens* its KD-tree query
  radius to $r_q = 0.05 + v_{\max}(\Delta t/2) \approx 0.0536$ AU with
  $v_{\max}=25$ km/s, which — under near-linear relative motion — provably brackets any
  encounter with $v_{\rm rel}\le 25$ km/s between two 12 h samples (it is within $r_q$
  at the nearer sample and enters the candidate set; refinement then finds the true
  minimum). The only encounters the cadence can miss have $v_{\rm rel}>25$ km/s, and a
  direct query of the catalogue's own velocity distribution shows **only 0.0012 %**
  exceed that (p99.99 = 20 km/s, max 68 km/s) — so the residual cadence-induced
  incompleteness is $\lesssim10^{-3}$ %, three orders of magnitude below the 0.70 %
  censoring term. The manuscript (§2.2, §3.3) has been corrected to describe the
  widened radius and state this bound. This turns an apparent hole into a bounded,
  near-zero budget term and *strengthens* the completeness thesis.
- **(Moderate) The reported distances are two-body except on the fragile subset.**
  Correctly disclosed, but the abstract/title say "real 3D close encounters"; a strict
  reader wants "candidate" nearer the front. The hybrid catalogue mitigates this for
  the fragile subset only.
- (Minor) Mutual gravitational deflection between the two encountering bodies is
  neglected in the geometry (both are test particles against the planetary+big-16
  field). Justified at these masses, but worth an explicit sentence.

**Verdict:** sound methods; the one major analytical concern (cadence) was closed this
pass. **Minor revision** on the remaining dynamics wording (the "candidate" framing and
mutual-deflection sentence below).

---

## C. Referee 2 — Statistics & the completeness claim (the novelty)

**Strengths.**
- The threshold-censoring analysis is genuinely careful: it recognises the one-sided
  "no false negatives" as *censoring* (a pair Kepler placed above 0.05 AU was never
  written), then measures the missing side directly by widening to 0.06 AU and
  re-refining the [0.05, 0.06) band under N-body — 0.70 % [95 % CI 0.59–0.83 %] cross
  downward. Framing it as an Eddington/Malmquist selection effect with a symmetric,
  scatter-dominated Δd (σ ≈ 0.37 mAU) is correct and well argued. This is publishable
  statistics.
- Reporting a **jackknife (external) σ** that is 2–3× the formal Fisher σ, and
  *demonstrating its necessity on the calibrators themselves* (Pallas/Vesta sit at
  |z|>3 under formal σ but |z|<3 under jackknife), is a strong, self-validating move.
- The identifiability criterion (measured / not_identifiable / bound via χ² curvature)
  is the honest response to the mass↔orbit degeneracy.

**Concerns.**
- **(Major) The prefilter-recall result (76 %) does not describe the delivered
  catalogue.** The paper itself states the pair-list prefilter was *skipped* at
  main-belt scale, so the |Δa| cut "did not actually drop these pairs from the frozen
  catalogue"; the 76 % measures damage the prefilter *would* cause "if applied at small
  N." Presenting a prospective what-if as one of three pillars of a *measured* budget
  for *this* catalogue is a framing problem. A referee would ask to either (i) demote
  it to a "recommendation for future runs" subsection, or (ii) keep it but state
  unambiguously in the abstract that it does not bound this freeze. As written the
  abstract lists "76 %" beside the two terms that *do* apply, inviting misreading.
- **(Moderate) The catalog-wide extrapolation to ~10⁵ censored encounters** is
  explicitly order-of-magnitude and rests on an N² scaling of the [0.05, 0.06) band.
  Acceptable if flagged (it is), but a referee may want the scaling justified or the
  claim softened to "10⁵–10⁶".
- (Minor) The 0.70 % is measured on 10,000 numbered bodies; a sentence on
  representativeness vs the full ~150k would pre-empt a question.

**Verdict:** the central novelty is real and defensible. After this pass the *measured
budget* has two measured-and-applicable terms (censoring, Kepler error) plus a
bounded-and-applicable cadence term (≲10⁻³ %), leaving only the prefilter-recall
framing to fix (measured but not applicable to this freeze). Tighten that one framing
point → **minor–moderate revision**, not reject: the science is sound.

---

## D. Referee 3 — Catalogue / data / reproducibility

**Strengths.**
- Provenance discipline is excellent: frozen MPCORB snapshot with SHA, catalogue SHA,
  configuration sidecar, explicit non-bit-reproducibility caveat. This is better than
  most catalogue papers.
- Three tiered products (Kepler candidate / hybrid / characterised) with per-row
  `refinement_method` on the hybrid is exactly what downstream users need.
- Every completeness number is tied to a named, rerunnable script.

**Concerns.**
- **(Major, blocking for publication) No DOI / archive yet.** A&A (and PSJ) require the
  data at CDS/VizieR and/or a Zenodo DOI before acceptance. Scaffolding now exists
  (`DATA_AVAILABILITY.md`, `zenodo_data_deposit.json`) but the deposit must be executed
  and the DOI inserted. Non-negotiable.
- **(Moderate) 15 GB across three parquet files** is awkward for CDS. The referee/
  editor will want a documented column schema (a ReadMe) and likely a reduced
  machine-readable table at CDS with the bulk on Zenodo.
- (Minor) Author list, affiliations, ORCID, acknowledgements, software-citation stack
  (astropy, rebound, ASSIST, polars, scipy) are placeholders. Standard but required.

**Verdict:** **accept conditional on the data deposit + author metadata**; the
engineering and provenance are above the bar.

---

## E. Comparison with analogous published work

*(Citations verified against publisher pages — A&A/IOP/Springer — in a dedicated
four-front literature search. ADS abstract pages render client-side and could not be
fetched, so volume/page/DOI come from the publishers directly.)*

**E.1 Encounter / close-approach catalogues.**

| Work | Citation | Scale | Completeness budget? |
|------|----------|------:|----------------------|
| Galád | A&A 370, 311 (2001) | 9,511 asteroids | qualitative only |
| Galád & Gray | A&A 391, 1115 (2002) | 24,599 asteroids | "not guaranteed all pairs found" — **unquantified** |
| Baer & Chesley | CeMDA 100, 27 (2008) | 21 masses; defines the two-body encounter-selection algorithm reused since | none |
| Baer, Chesley & Matson | AJ 141, 143 (2011) | 26 masses | "several successful events… missed" — **unquantified** |
| Fuentes-Muñoz et al. | AJ 170, 353 (2025) | **975,673** test asteroids screened (dCA < 0.05 AU) | none |
| Bernstein, Najafi & Gomes | arXiv:2504.01184 (LSST, 2025) | 6,839 encounters | "not exhaustive" — qualitative |

**Verified bottom line:** no prior published asteroid–asteroid mutual-encounter list
reports a *measured* completeness budget. Several authors *know* their lists are
incomplete and say so **in words** (Galád & Gray 2002; Baer et al. 2011; Bernstein et
al. 2025) — none converts it into numbers. **The novelty claim holds**, with two
framing refinements the paper should adopt: (i) phrase it as "first *for
asteroid–asteroid mutual-encounter catalogues*," since completeness budgeting is
routine in adjacent domains (impact-monitoring; Gaia stellar-encounter rates); (ii) the
raw scale (72.2 M) independently dwarfs the nearest prior effort (FM 2025's ~976 k
screened test asteroids). Note: the "21 asteroids" work is **Baer & Chesley 2008, not
Goffin** — the manuscript's Goffin citation is for a *simultaneous* determination
(A&A 565, A56: 132 acceptable masses, 49 new), which is correct.

**E.2 Mass determination — state of the art.**
Baer & Chesley (2008/2011/2017: 21/26/30 masses), Goffin (2014: 132 masses), Konopliv
et al. (2006/2011: the reference masses for the biggest perturbers via Mars ranging),
Michalak (2000/2001), and the direct competitor **Fuentes-Muñoz et al. 2025** (AJ 170,
353): **232 masses at SNR > 3, 159 at > 5, 77 at > 10**, same mechanism and same Gaia
FPR data. On raw count the manuscript's "16 perturbers / 10 identifiable" is
**confirmatory, not catalogue-competitive** — the paper says so. But the two
methodological novelties are **confirmed genuine against the competitor**:
- **Full intra-transit / FOV-block along-scan covariance**: FM 2024/2025 use only a
  *scalar* per-transit systematic floor `C = C_rnd + C_sys/N`; Baer et al. 2011/2017
  model same-observatory correlations for ground data but not Gaia's FOV blocks. **No
  deflection-mass paper models the full block covariance** the manuscript uses
  (ICC = 0.32, √1.66 inflation).
- **External (jackknife/bootstrap) error bars**: every deflection-mass paper uses
  *scaled formal covariance* (Baer & Chesley: variational re-runs; FM: scaled formal;
  the exception is Siltala & Granvik 2020, A&A 633, A46, Bayesian MCMC). Crucially,
  **Siltala & Granvik explicitly argue LSQ literature uncertainties are underestimated /
  non-Gaussian** — that standing critique is *exactly* what the manuscript's jackknife σ
  answers. This is the sharpest way to sell the error methodology to a referee.

**E.3 Gaia context** (Tanga et al. 2023, A&A 674, A12; FPR orbits David et al. 2023,
A&A 680, A37): ~156 k asteroids ("~150 k" is a fair round-down); **along-scan precision
sub-mas (best ~0.25 mas), across-scan ~1 arcsec** — a ~3-order AL/AC asymmetry that is
the published physical motivation for the anisotropic along-scan covariance. **Window
correction applied to the manuscript:** the DR3 SSO nominal-ops window is 2014 Aug 5 –
2017 May 28; the propagation window (2014 Jul 25 – 2017 May 28) brackets it (§2.1 now
says so). The mass layer's FPR astrometry spans a longer arc (to 2020 Jan) — the paper
correctly says "full arc of Gaia FPR astrometry," so no conflation there.

> The comparison confirms the **shape** of the contribution: novel as a *dataset +
> methods* paper (first encounter catalogue with a measured budget; the first
> deflection-mass treatment with FOV-block covariance + external jackknife σ, which
> directly answers the Siltala & Granvik critique), not as a *discovery* paper (no new
> mass, no new large–large encounter — stated honestly).

---

## F. Editor's synthesis and verdict

**What kind of paper is this?** A **dataset + methods** paper. Judged as such it is
solid, honest, and technically careful, with a genuinely novel framing (a measured
completeness budget for an asteroid–asteroid encounter catalogue) that the literature
comparison supports. Judged as a **discovery** paper it would fail — and the authors
correctly do not frame it that way.

**Is the central claim ("measured completeness budget") delivered?** *After this pass,
essentially yes.*
- Kepler-vs-N-body error: **measured and applicable** ✔
- Threshold censoring (0.70 %): **measured and applicable** ✔
- Coarse-cadence miss rate: **bounded ≲10⁻³ %** this pass (widened-radius bracketing +
  the 0.0012 % v_rel>25 km/s tail); §2.2/§3.3 corrected ✔
- Prefilter recall (76 %): measured but **not applicable to this freeze** (prefilter
  skipped) — the remaining framing item; reframe as prospective.

**Consolidated verdict: MINOR–MODERATE REVISION** (clear path to acceptance at A&A or
PSJ). Remaining before acceptance:

1. ~~Measure the coarse-cadence loss~~ — **done this pass**: widened-radius bracketing
   guarantees capture to 25 km/s; residual ≲10⁻³ %. (Was the decisive item.)
2. **Reframe the 76 % prefilter recall** as a prospective recommendation, and make the
   abstract distinguish terms that bound *this* catalogue (Kepler error, censoring,
   cadence) from the term that bounds *future prefiltered* runs (prefilter recall).
   (Referee 2 — the main remaining analytical wording item.)
3. **Execute the data deposit** — and treat it as the **single highest-probability
   referee/editor trigger**. A&A's data policy is *mandatory and contractually
   enforced*: a 72 M-row catalogue must go to CDS/VizieR as documented machine-readable
   tables (compliant ReadMe, units), typically with a reduced table at CDS + the bulk on
   Zenodo. CDS independently checks consistency and can hold the paper at proofs.
   Scaffolding is in `docs/paper/` (`DATA_AVAILABILITY.md`, `zenodo_data_deposit.json`);
   fill author metadata + software citations. (Referee 3, blocking-but-mechanical.)
4. **Anticipate an injection–recovery demand on the completeness budget.** The
   censoring measurement (re-refining the [0.05, 0.06) band) and the FM-2025 literature
   recovery (§3.4) are recovery experiments, but an A&A referee for a completeness-claim
   paper will likely ask for a *synthetic injection–recovery* (inject known encounters,
   report recovered fraction and both false-negative **and false-positive** rates). The
   pieces largely exist; framing them as a single injection–recovery section would
   pre-empt the most probable *scientific* revision request. (New, from the literature
   search of A&A norms.)
5. **Frame the mass-error methodology against Siltala & Granvik (2020).** State
   explicitly that jackknife externals answer the standing critique that LSQ literature
   σ are underestimated/non-Gaussian, and that the FOV-block covariance is absent even
   from the direct competitor (FM 2025). This converts "confirmatory" into a defensible
   *methods* advance. (§E.2.)
6. Minor (partly done this pass): "candidate" now in the abstract; mutual-deflection and
   10k-body representativeness sentences added (§2.3, §3.2); still to do — soften/justify
   the 10⁵ extrapolation. DR3 window start clarified (§2.1); "231" → "over 230 (232 at
   SNR>3)".

**What would make it a strong accept:** items 1–3 done. The methods (FOV-block
covariance, jackknife σ, identifiability, measured censoring) and the provenance
discipline are already above the typical bar; the catalogue is a real community
resource. The honesty about null results (no new mass/encounter) is a credit, not a
liability, for a dataset paper — provided the abstract sells the *budget + methods* as
the result, which it now largely does.

**Probable outcome if resubmitted with 1–4 addressed:** accept at A&A as a
catalogue/methods paper.

---

## G. Actionable checklist distilled for the authors

- [x] Quantify the coarse-cadence miss rate — **done this pass**: the widened query
      radius brackets v_rel ≤ 25 km/s and only 0.0012 % of catalogued encounters exceed
      that, so cadence loss is ≲10⁻³ %. §2.2/§3.3 corrected.
- [ ] Reframe §3.3 prefilter recall as prospective; fix the abstract to separate
      "bounds this freeze" (censoring, Kepler error) from "bounds future runs"
      (prefilter). 
- [ ] Execute Zenodo/CDS deposit; insert DOI (scaffolding in `docs/paper/`).
- [ ] Author metadata, acknowledgements, software citations.
- [ ] Minor wording: "candidate" framing; soften 10⁵ extrapolation; mutual-deflection
      and representativeness sentences.
