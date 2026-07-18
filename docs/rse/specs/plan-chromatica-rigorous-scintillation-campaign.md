# Implementation Plan: Chromatica rigorous scintillation campaign

**Date:** 2026-07-17
**Author:** AI Assistant
**Status:** Complete — implemented, validated, pushed, and documented in draft PR #200
**Related Documents:**
- [Research: rigorous scintillation campaign](research-chromatica-rigorous-scintillation-campaign.md)
- [Experiment: fit and uncertainty architecture](experiment-chromatica-scintillation-fit-uncertainty.md)

## Overview

Replace the provisional DSA point-loading path in PR #200 with a reproducible,
common CHIME/FRB and DSA-110 measurement contract. Telescope-specific code may
prepare the dynamic spectra, masks, windows, and harmonic exclusions, but both bands
will use the same weighted Lorentzian selector, bootstrap semantics, component
reporting, and fail-closed qualification.

**Goal:** Produce requalified bandwidth and modulation records, then regenerate the
cross-band power-law result only from points that pass every required gate.

## Current State Analysis

- `scintillation/scint_analysis/revalidation.py:277-425` provides the weighted
  one-to-three Lorentzian central estimator.
- `scintillation/scint_analysis/window_refit.py:112-390` implements the accepted
  CHIME-specific fit and window campaign but does not expose ACF arrays to a common
  postprocessor.
- `analysis/scintillation-dsa-lorentzian-2026-07-07/run_dsa_lorentzian_fits.py:1338-1555`
  computes DSA controls, but its final non-CHIME status does not fail closed.
- `analysis/chromatica-cross-band-scintillation-2026-07-17/fit_cross_band.py`
  currently accepts the old campaign schema and uses DSA conditional covariance
  errors without a DSA systematic.

## Desired End State

- A pure, tested common fitter reports ordered components, `m_narrow`, `m_broad`,
  and `m_total`, block-bootstrap intervals, model-count stability, alternative-shape
  sensitivity, and explicit admission reasons.
- A Chromatica runner evaluates central and variant ACFs for both telescopes,
  including on/off windows, first retained lag, fit span, and 2/3/4 equal-S/N
  partitions.
- Matched-data injection uses the real subband channelization, mask, off-pulse noise,
  and burst-window duration.
- Every central subband is individually fail-closed on normalization, fit quality,
  off-pulse null, low-lag stability, bootstrap/model stability, variant stability,
  alternative shape, and injection recovery.
- Cross-band fitting consumes only the rigorous schema. If fewer DSA points qualify,
  the absence of a joint fit is recorded as a gate result rather than bypassed.

## What We're NOT Doing

- Rewriting the telescope-specific calibration or RFI masking algorithms.
- Promoting `m_total` when a broad component is unresolved or inconsistently selected.
- Updating the `Faber2026` parent submodule pointer in this task.
- Merging PR #200 before the scientific outputs and figures are reviewed.

## Implementation Approach

1. **Common central estimator:** retain the weighted BIC plus corroborating F-test
   selector. The architecture experiment found no benefit in replacing it with the
   unweighted CHIME optimizer.
2. **Layered uncertainty:** combine block-bootstrap half-interval, central-fit
   covariance, and variant half-range in quadrature; retain each term separately.
3. **Component semantics:** sort by width; the narrowest component is the bandwidth
   candidate. `m_broad` is reported only for a stable resolved second component, and
   `m_total = sqrt(sum(m_i^2))` is eligible only when all components are eligible.
4. **Matched-data injection:** add a known thin-screen source spectrum to real
   off-pulse spectra of the same duration, preserving channel masks and the production
   ACF path. Recovery and coverage gates use fixed seeds.
5. **Fail closed:** `True` is required. Missing or inconclusive gates are failures,
   and every failure reason is serialized.

## Implementation Phases

### Phase 1: Common fit and report contract

- [x] Add failing unit tests for component ordering, modulation semantics, correlated
  block resampling, deterministic seeds, and fail-closed missing gates.
- [x] Run the focused tests and confirm they fail because the module does not exist.
- [x] Implement `scintillation/scint_analysis/rigorous_campaign.py` with the weighted
  central fit, generalized-Lorentzian sensitivity, moving-block bootstrap, uncertainty
  combination, and qualification functions.
- [x] Run the focused tests and confirm they pass.

### Phase 2: Telescope adapters and campaign variants

- [x] Add failing adapter tests for ACF export and rigorous-schema construction.
- [x] Expose the required CHIME ACF arrays from `window_refit.refit` without changing
  its existing result fields.
- [x] Implement
  `analysis/chromatica-cross-band-scintillation-2026-07-17/run_rigorous_campaign.py`
  to prepare central/variant ACFs for CHIME and DSA and run the shared contract.
- [x] Include 2/3/4 subband partitions as campaign-level consistency controls; use
  fixed four-subband boundaries for pointwise window and fit-policy systematics.

### Phase 3: Matched-data injection and scientific products

- [x] Add analytic tests showing recovery of a known Lorentzian width and modulation
  in noise-free data within discretization tolerance.
- [x] Add matched off-pulse injection/recovery for every central subband using real
  channel masks, off-pulse slices, channel widths, and burst duration.
- [x] Run the fixed-seed campaign; serialize full provenance, thresholds, commands,
  hashes, package versions, and all failed gates.
- [x] Render per-subband gate/fit figures and a modulation-eligibility figure without
  plot titles.

### Phase 4: Fail-closed cross-band regeneration

- [x] Add failing tests that reject the legacy DSA schema and refuse a joint result
  when fewer than two qualified DSA points remain.
- [x] Update `fit_cross_band.py` to consume the rigorous result, propagate the full
  width uncertainty, and record when the joint law is unavailable.
- [x] Regenerate JSON, CSV, PNG, and PDF results; update the analysis README and figure
  review record with the new scientific status.

### Phase 5: Validation and publication

- [x] Run focused and existing regression tests, Ruff on touched Python, and the
  campaign from the locked environment.
- [x] Repeat the result-generating commands in an isolated clean environment and
  compare deterministic hashes or declared numeric tolerances.
- [x] Inspect every new figure at full resolution and record the review.
- [x] Write implementation and validation reports, run `agent-closeout-check`, commit
  only task-scoped paths, push, and update draft PR #200.

## Success Criteria

### Automated Verification

- The new focused tests pass.
- Existing window-campaign, ACF-covariance, revalidation, and cross-band tests pass.
- A fixed-seed rerun reproduces admitted points and intervals within `1e-8` relative
  tolerance.
- No legacy DSA point can be silently accepted by the cross-band loader.
- Every output point contains explicit gate booleans, failure reasons, uncertainty
  components, input hashes, and seed provenance.

### Scientific Correctness

- A synthetic Lorentzian ACF recovers injected gamma and modulation within 5% in the
  noise-free resolved regime.
- Matched-data injection median relative gamma bias is at most 15%, empirical 68%
  interval coverage is at least 60% for the finite trial count, and recovery succeeds
  in at least 80% of trials.
- Bootstrap preferred-model fraction is at least 70% and variant preferred-model
  fraction at least 60% for an admitted modulation component. A narrow bandwidth may
  remain stable across a model-count change, but modulation remains conditional on the
  decomposition and is excluded.
- Alternative-shape and fit-policy widths differ from the central width by no more
  than 35% for an admitted point.
- Physical modulation requires `0 < m <= 1.2`; `m_total` additionally requires every
  selected component to pass its component gate.

## Risk Assessment

1. **Insufficient off-pulse slices:** return an inconclusive injection/null gate and
   exclude the point; never synthesize a pass.
2. **Subband boundaries move across partitions:** treat 2/3/4 partitions as a
   campaign-level consistency diagnostic, not a one-to-one point uncertainty.
3. **Runtime:** cache prepared ACF records and separate cheap fit-policy variants from
   expensive window/subband preprocessing; keep fixed trial counts in provenance.
4. **No qualified points:** publish a scientifically useful failed-qualification
   result and any independently qualified single-band characterization; do not force a
   joint or single-band slope when its support gate fails.

## Open Questions

None. The user approved autonomous full implementation on 2026-07-17.

## Review History

### Version 1.0 — 2026-07-17
- Approved direct-mode implementation plan derived from the code audit and fit
  architecture experiment.
