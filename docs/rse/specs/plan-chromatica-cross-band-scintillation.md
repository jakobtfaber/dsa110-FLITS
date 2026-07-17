# Plan: Chromatica cross-band scintillation fit

Date: 2026-07-17
Research basis: `docs/rse/specs/research-chromatica-cross-band-scintillation.md`

## Goal

Produce a reproducible CHIME+DSA-110 power-law characterization for Chromatica
without admitting the DSA subband that fails the off-pulse instrumental null.

## Phase 1: Lock the fresh DSA measurement

- Run the existing fresh Lorentzian driver for Chromatica from the staged DSA `.npz`.
- Preserve the selected result JSON, component CSV, report, and diagnostic figures.
- Verify that the selected subband candidate and per-subband artifact-control fields
  are present.

Success: the run completes without fallback to YAML `stored_fits`, and its input and
output provenance is retained under the analysis directory.

## Phase 2: Implement the cross-band fitter

- Parse accepted CHIME narrow-component measurements and their three uncertainty
  terms.
- Parse the selected DSA `gamma_1` components and enforce per-subband quality,
  off-pulse-null, and low-lag-stability gates.
- Fit the power law in log space using measurement errors.
- Fit a second likelihood with a non-negative intrinsic log-scatter term.
- Calculate chi-square, degrees of freedom, p-value, covariance-based uncertainties,
  residuals, and leave-one-DSA-out sensitivities.
- Emit JSON, CSV, and title-free PDF/PNG figures using mathematical axis notation.

Success: invalid/non-positive inputs fail closed; the instrumental 1459.62 MHz DSA
point is recorded as excluded; outputs identify the primary model and why it was
chosen.

## Phase 3: Verification and publication

- Unit-test synthetic parameter recovery, gate behavior, and malformed-input failure.
- Re-run the analysis from the locked inputs and compare deterministic numerical
  outputs.
- Inspect the diagnostic and cross-band figures visually.
- Run the repository closeout gate with the dirty-state/restart packet.
- Commit only this analysis, its research/plan documents, and generated products;
  preserve the parent-repo submodule pin and unrelated recovered artifacts.
- Push the focused branch and open a draft pull request for scientific review.

Success: tests and closeout pass, the task-scoped commit is published, and the PR
states that the exact no-scatter law is rejected even though a cross-band slope with
extra scatter is measurable.
