# Validation Complete

> Validated against `docs/rse/specs/plan-chromatica-cross-band-scintillation.md` /
> `docs/rse/specs/implement-chromatica-cross-band-scintillation.md` at commit
> `ae67bdf` on 2026-07-17.

## Overall Status: Ready for a draft scientific-review PR

## Summary

- Phases: phases 1 and 2 fully implemented; phase 3 verification complete and
  publication pending this validation.
- Focused automated checks: 4 passing check groups, 0 task-scoped failures.
- Broader driver checks: 14 passing, 1 unrelated pre-existing roster-artifact failure.
- Manual testing: both relevant figures visually inspected; scientific promotion of
  the 1396 MHz DSA point remains a human review item.
- Critical issues: 0.
- Important issues: 1 scientific follow-up (DSA subband scatter/systematics).

## Implementation Status

### Phase 1: Lock the fresh DSA measurement

**Status:** Fully implemented

- Fresh ACF/Lorentzian run from `chromatica.npz`: complete.
- Selected four-subband result, aggregate result, CSV/report, and figures preserved:
  complete.
- Candidate-selection and per-subband artifact-control fields inspected: complete.

### Phase 2: Implement the cross-band fitter

**Status:** Fully implemented

- CHIME admission and three-term uncertainty propagation: complete.
- DSA `gamma_1` admission with component, off-pulse, and low-lag gates: complete.
- Formal and intrinsic-scatter likelihoods: complete.
- Goodness of fit, covariance errors, residuals, and leave-one-DSA-out sensitivity:
  complete.
- JSON/CSV/PDF/PNG outputs and title-free project figure styling: complete.
- Near-zero-scatter boundary uncertainty handled without infinity/overflow: complete.

### Phase 3: Verification and publication

**Status:** Verification complete; publication follows this report

- Unit tests, lint, format, deterministic JSON rerun, and visual review: complete.
- Closeout, focused commit, push, and draft PR: pending at validation time.

## Automated Verification Results

### Passing checks

- `pytest .../test_fit_cross_band.py` — 5 passed.
- `pytest .../test_fit_cross_band.py .../test_driver_guards.py -k 'not pbf_loader'`
  — 14 passed, 1 deliberately deselected.
- `ruff check` — all checks passed.
- `ruff format --check` — both files already formatted.
- Two consecutive fitter runs — identical result JSON SHA-256
  `bf97c437f8dbff0cb6d9108c7f635f8bcae2b2b9256a7aa18c975e721dba6439`.
- Fresh numerical result — alpha `3.632970 +/- 0.392372`, gamma at 1 GHz
  `0.328985 MHz`, intrinsic log scatter `0.323774`; formal fit chi-square/dof
  `68.196/5`.

### Broader check with an unrelated failure

- Full `test_driver_guards.py` run — 14 passed and 1 failed. The failure is
  `test_pbf_loader_follows_locked_roster_and_excludes_gate_failures`: the locked
  Mahi roster points to `mahi_joint_fit_C1D2.json`, which is absent in this recovered
  commit. The driver correctly fails closed and returns `None`.
- Root cause: missing pre-existing Mahi scattering/PBF roster artifact, unrelated to
  Chromatica and to the files changed here.
- Regression assessment: no regression. `git diff` confirms this task does not modify
  the DSA driver or its tests, and the Chromatica run does not require the missing Mahi
  overlay.
- Recommendation: repair the Mahi roster/artifact pairing in a separate task rather
  than weakening the fail-closed loader in this fit.

## Code Review Findings

### What matches the plan

- Every accepted/rejected point and the reason for rejection is retained.
- The 1459.62 MHz DSA point is excluded specifically because its off-pulse null fails.
- The exact power law is not promoted: its p-value is `2.43e-13`.
- The primary result includes intrinsic log scatter and exposes sensitivity to each
  DSA point.
- Input JSON, raw DSA data, and fitter hashes are recorded.

### Deviations from the plan

- None in task scope.

### Potential issues

- The current DSA driver supplies Lorentzian covariance errors but not a DSA window-
  campaign systematic analogous to CHIME's. The primary intrinsic-scatter model
  absorbs observed excess dispersion statistically, but this is not a substitute for
  a dedicated DSA systematic campaign.
- The 1395.89 MHz point drives much of the excess scatter: leaving it out changes alpha
  from `3.63 +/- 0.39` to `3.12 +/- 0.19` and drives fitted intrinsic scatter to the
  zero boundary. This sensitivity is disclosed in the result JSON.

## Manual Testing Required

1. **Scientific promotion decision**
   - Review the DSA ACF panel and the cross-band fit, focusing on the accepted
     1395.89 MHz component and the structured 1321 MHz residuals.
   - Decide whether the result should remain a diagnostic/campaign result or enter the
     manuscript. No manuscript claim is changed by this PR.

## Recommendations

### Critical (must fix before draft PR)

- None.

### Important (should fix before manuscript promotion)

- Run a DSA window/subband systematic campaign and quantify selection variability,
  especially for the 1395.89 MHz component.

### Nice to have

- Add a profile-likelihood or bootstrap interval once the DSA systematic model is
  defined.

### Follow-up work

- Repair the unrelated Mahi PBF roster artifact expected by the broader driver test.

## References

- Plan: `docs/rse/specs/plan-chromatica-cross-band-scintillation.md`
- Implementation: `docs/rse/specs/implement-chromatica-cross-band-scintillation.md`
- Research: `docs/rse/specs/research-chromatica-cross-band-scintillation.md`
