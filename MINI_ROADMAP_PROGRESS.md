# Mini-Roadmap Progress: Davida Mass Fit

> Tracking file for the 5-day proof-of-concept from ROADMAP_TO_MASSES.md.
> Update task status as work progresses so the next session can resume cleanly.
>
> **Goal**: Recover the mass of (511) Davida from its 2014-11-19 encounter with
> asteroid 2003_sm90 using Gaia DR3 astrometry. Compare against Goffin (2014)
> value of ~3.5e19 kg.
>
> **Branch**: `feat/mini-roadmap-davida-mass`

---

## Status legend
- 🔲 not started
- 🟡 in progress
- ✅ done
- ❌ failed / blocked
- ⏭️ skipped (with rationale)

---

## Day 1-2: Coordinate transforms + N-body wrapper

### T1.1 — Module skeleton ✅
- [x] `src/astrometry/__init__.py`
- [x] `src/astrometry/transforms.py`
- [x] `tests/test_astrometry_transforms.py`

### T1.2 — Ecliptic → equatorial ✅
- [x] `ecliptic_to_equatorial(xyz_ecl)` — J2000 obliquity rotation
- [x] Unit test: X-axis invariant, Y-axis rotates correctly
- [x] Inverse: `equatorial_to_ecliptic` with roundtrip test

### T1.3 — Heliocentric → barycentric ✅
- [x] `sun_barycentric_au(jd_tdb)` — uses astropy `get_body_barycentric`
- [x] `heliocentric_to_barycentric_icrs(pos_helio, jd_tdb)` — adds Sun's barycentric position
- [x] Test: Sun's barycentric position has |r_Sun| < 0.01 AU

### T1.4 — Light-time correction ✅
- [x] `light_time_iterate(target_pos_func, jd_tdb_obs, gaia_xyz, max_iter=3, tol_seconds=1.0)`
- [x] Iterative, converges in 2-3 steps
- [x] Tested with stationary target (τ = d/c) and moving target

### T1.5 — Stellar aberration ✅
- [x] `stellar_aberration(line_of_sight, observer_velocity)` — first-order
- [x] Test: zero velocity → no aberration
- [x] Test: 30 km/s perpendicular → ~20 arcsec deflection
- [ ] **TODO**: still need to compute Gaia's velocity vector from x_gaia(t) time series.
      For now the function expects observer_vel as input. Will be wired in T2.1.

### T1.6 — Sanity test: predict RA/Dec for Ceres vs Horizons 🟡
- [x] `scripts/sanity_check_transforms.py` — full chain (Kepler + transforms)
- [x] `scripts/sanity_check_transforms_only.py` — transforms only (Horizons vectors as input)
- [x] **Result with full chain (Kepler 2-body from 2012 epoch)**: 17 arcsec residual.
      Dominated by Kepler drift from missing planet perturbations, NOT transform bugs.
- [x] **Result transforms-only (Horizons vectors)**:
      - Without light-time: 11 arcsec  (light-time matters!)
      - With light-time:    1 arcsec   (✓ light-time validated)
      - With light-time + aberration:  12 arcsec  (aberration makes it worse)
- [x] **Discovery**: Gaia DR3 SSO reports (ra, dec) in barycentric astrometric ICRS,
      with stellar aberration already removed by the Gaia pipeline.
      So we must NOT apply our own aberration on top.
- [x] Final operating residual: **~1 arcsec** (acceptable for step detection;
      remaining error likely from Horizons-vector spline interpolation +
      solar gravitational deflection not yet modelled). Good enough for step
      detection; may need refinement for precision mass fit.

### T1.6b — Open question / TODO 🔲
- [ ] Investigate the residual 1 arcsec: is it interpolation, solar deflection,
      something else? Probably worth re-checking after T1.7 once we use N-body
      positions instead of Horizons-vector splines.

### T1.7 — N-body wrapper with configurable perturber 🔲
- [ ] `src/propagate/nbody_perturber.py`
- [ ] `propagate_with_perturber(target_elements, perturber_number, perturber_mass_kg, t_start_jd, t_end_jd, step_days, include_planets=True, include_big4=True) → (jd_array, xyz_array)`
- [ ] Wraps REBOUND; adds the perturber as a 4th additional point mass (beyond big-4)
- [ ] Loads perturber's orbital elements from MPCORB
- [ ] Test: with mass=0, output ≈ baseline (no perturber); with mass=known, recover the deflection

### T1.8 — End-to-end test: Davida perturbing test particle 🔲
- [ ] Set up: (511) Davida as perturber with mass = 3.5e19 kg (Goffin)
- [ ] Target: a test particle near 2003_sm90's orbit
- [ ] Propagate ±90 days through 2014-11-19 encounter
- [ ] Verify: total angular shift on the test particle matches expected ~3 mas

---

## Day 3: Forward model + fit machinery

### T2.1 — Forward model 🔲
- [ ] `forward_model(elements, perturber_mass, perturber_number, epochs_jd_tdb, gaia_xyz_array) → (ra_deg, dec_deg)`
- [ ] Internally: propagate → light-time iterate → aberration → ICRS RA/Dec
- [ ] Vectorised over epochs

### T2.2 — Residual function for scipy 🔲
- [ ] `residuals(params, args) → 1D array of (Δra·cos(dec), Δdec) in mas`
- [ ] params = (a, e, i, Omega, omega, M0, log10_mass)
- [ ] Bounds: mass > 0 (work in log space), eccentricity ∈ [0, 1), inclination ∈ [0, π)

### T2.3 — Initial guess 🔲
- [ ] Orbital elements: from MPCORB
- [ ] Mass: from diameter via ρ=1.5 g/cm³

### T2.4 — Fit script CLI 🔲
- [ ] `scripts/fit_perturber_mass.py --perturber N --target M --date YYYY-MM-DD`
- [ ] Output: fitted mass + sigma, chi2_red, residual CSV
- [ ] Method: scipy.optimize.least_squares (trf with bounds, jacobian-based)
- [ ] Uncertainty: Σ = (Jᵀ J)⁻¹ · χ²_red (1σ via diagonal)

---

## Day 4: Apply to Davida

### T3.1 — Run on (511) Davida + 2003_sm90 🔲
- [ ] Encounter date: 2014-11-19
- [ ] Window: ±180 days
- [ ] Expected mass: ~3.5e19 kg (Goffin 2014)

### T3.2 — Save diagnostics 🔲
- [ ] `data/output/davida_fit.csv` — residuals per epoch
- [ ] `data/output/davida_fit_summary.json` — fitted params + uncertainties

### T3.3 — Inspect 🔲
- [ ] Check χ²_red — should be ~1 if errors well-estimated, ~few if not
- [ ] Check post-fit residuals — should be at mas level
- [ ] Check that mass posterior is constrained (not unbounded)

---

## Day 5: Validation

### T4.1 — Compare against Goffin 2014 🔲
- [ ] Pull Goffin's published mass for Davida (~3.5e19 kg)
- [ ] Compare fitted ± σ vs published
- [ ] If within factor 2 → mini-roadmap succeeded ✅
- [ ] If off by factor >5 → debug; method has fundamental issue

### T4.2 — Document result 🔲
- [ ] Write `MINI_ROADMAP_RESULT.md` with:
  - Fitted mass and uncertainty
  - χ²_red, n_obs, n_iterations
  - Comparison to Goffin
  - List of caveats discovered during implementation
  - Decision: proceed to full catalog OR fix bugs first

### T4.3 — Open PR 🔲
- [ ] Branch `feat/mini-roadmap-davida-mass` → main
- [ ] PR description summarises result and next steps

---

## Notes / decisions log

(Add notes here as we go — bugs found, design decisions, unexpected things)

- 2026-05-19: Mini-roadmap started. Branch created. Beginning with T1.1.

---

## Quick reference for resuming

If a future session picks this up:

1. `git checkout feat/mini-roadmap-davida-mass`
2. Read this file top-to-bottom to see what's done
3. Continue from the first 🔲 task
4. Update status as you go (🔲 → 🟡 → ✅)
5. When all tasks done, generate `MINI_ROADMAP_RESULT.md` and open PR

Key files that already exist (don't re-implement):
- `src/propagate/nbody.py` — REBOUND wrapper (already supports big-4)
- `src/propagate/kepler.py` — Kepler 2-body
- `src/ingest/mpcorb.py` — orbital element loader
- `scripts/check_gaia_observations.py` — TAP query pattern (epoch in days since J2010 TCB!)
- `scripts/demo_ate_clean.py` — Horizons ephemerides as observer (good reference pattern)

Key facts to remember:
- Gaia `sso_observation.epoch` = days since J2010.0 TCB (NOT a JD)
- Gaia observer location code in Horizons: `500@-139479`
- TAP sync limits to 2000 rows — use `launch_job_async`
- Davida MPCORB number: 511; 2003_sm90 number: 115180
