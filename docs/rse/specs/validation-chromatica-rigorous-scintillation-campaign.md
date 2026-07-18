# Validation Complete

> Validated against `plan-chromatica-rigorous-scintillation-campaign.md` /
> `implement-chromatica-rigorous-scintillation-campaign.md` at commit `169decb`
> on 2026-07-17.

## Overall Status: Ready

## Summary

- Phases: 5 of 5 implementation/validation phases complete; publication closeout is
  the only remaining operational step.
- Automated checks: 7 passing commands, 0 unresolved failures.
- Manual testing: 4 figures reviewed at original resolution; 0 items remain.
- Critical issues: 0.
- Important issues: 0.
- Scientific admission: 0/4 CHIME/FRB and 0/4 DSA-110 bandwidths; no modulation
  indices and no cross-band exponent are qualified.

## Implementation Status

### Phase 1: Common fit and report contract

**Status:** Fully implemented.

- The common weighted selector, component ordering, separate width/modulation
  eligibility, explicit-True gates, block bootstrap, named uncertainty terms, and
  generalized narrow-shape sensitivity are present in
  `scintillation/scint_analysis/rigorous_campaign.py`.
- Eight focused tests cover deterministic resampling, fail-closed missing gates,
  component semantics, known-truth recovery, two-scale alternative-shape behavior,
  and uncertainty combination.

### Phase 2: Telescope adapters and campaign variants

**Status:** Fully implemented.

- `window_refit.refit` retains its existing fields and additively exposes JSON-ready
  ACF data plus opt-in prepared objects.
- The analysis runner evaluates central and fixed-boundary window/fit-policy variants
  for both telescopes and serializes 2/3/4 equal-S/N partitions as campaign-level
  diagnostics.

### Phase 3: Matched-data injection and products

**Status:** Fully implemented.

- Every central subband is injected into real off-pulse spectra with its real channel
  slice, mask, channel width, source normalization, and burst duration.
- Seed, trials, recovery bias, final-interval coverage, and all thresholds are recorded.
- Three rigorous diagnostic figures were generated, hashed, and reviewed.

### Phase 4: Fail-closed cross-band regeneration

**Status:** Fully implemented.

- The production loader rejects legacy per-telescope JSON.
- A joint result requires at least two admitted points in each telescope.
- The regenerated result explicitly reports the joint fit unavailable; no slope or
  intrinsic-scatter envelope is produced.

### Phase 5: Validation

**Status:** Fully implemented.

- Focused/regression, optional-evidence, Ruff, diff-integrity, isolated reproducibility,
  and visual-review checks all completed.

## Automated Verification Results

### Passing checks

- `uv run --frozen pytest -q ...` on the common campaign, window campaign,
  multicomponent selector, artifact guards, revalidation, and cross-band tests:
  **66 passed in 7.84 s**.
- `uv run --frozen --extra nested pytest -q
  scintillation/scint_analysis/tests/test_acf_covariance.py`:
  **3 passed in 7.22 s**.
- `uv run --frozen ruff check` on the new/common campaign, covariance wrapper,
  campaign runner, cross-band fitter, and their tests: **all checks passed**.
- `git diff --check`: **passed**.
- `agent-closeout-check --packet /tmp/flits-rigorous-closeout.json`: **passed**;
  all dirty paths task-scoped, no runtime restart required.
- Locked-environment campaign run: **0/4 CHIME/FRB, 0/4 DSA-110 admitted**, command
  exited zero and wrote the full fail-closed record.
- `uv run --isolated --frozen` campaign in `/tmp/flits-rigorous-clean.DPMnSM`:
  `campaign_json_exact_match=true` and `figure_png_exact_match=true`.

### Investigated non-code failure

The first broad regression invocation ran `test_acf_covariance.py` without its declared
`nested` extra and produced one `ModuleNotFoundError: dynesty`; the other 60 tests in
that invocation passed. `pyproject.toml` declares `dynesty` under the `nested` extra.
Rerunning the complete file with `--extra nested` passed all three tests. This was an
environment invocation error, not a product failure, and is resolved in the final
verification command above.

## Code Review Findings

### What matches the plan

- The ACF calculator and telescope-specific preprocessing remain shared/reused rather
  than duplicated.
- Central fit choice follows the completed architecture experiment.
- Bootstrap draws reselect the component model and preserve contiguous lag residuals.
- Off-pulse and low-lag controls are per-subband and require explicit pass verdicts.
- Matched injections use observed noise/masks rather than ideal diagonal Gaussian ACF
  noise.
- `m_narrow`, `m_broad`, and `m_total` have independent eligibility and uncertainty
  fields; no conditional value is silently promoted.
- The cross-band result can no longer consume the legacy DSA result by accident.
- Result and figure provenance are sufficient to reproduce the exact artifacts.

### Deviations from the initial plan

1. **Width versus model-count stability were separated.**
   - Reason: a narrow width can be stable when a broad component appears/disappears,
     whereas modulation is conditional on the decomposition.
   - Impact: bandwidth qualification uses width stability; modulation additionally
     requires bootstrap and variant model-count stability.
   - Assessment: acceptable and more physically faithful.

2. **Injection coverage uses the final width interval.**
   - Reason: conditional covariance was already demonstrated to under-cover and is
     only one term of the declared uncertainty.
   - Impact: the fixed 60% coverage threshold tests the interval actually reported.
   - Assessment: required correction; the conditional-covariance version would test
     the wrong object.

3. **The campaign produced no single-band characterization.**
   - Reason: zero points in both bands passed every width gate.
   - Impact: the code reports no CHIME-only or cross-band law instead of forcing a fit.
   - Assessment: correct fail-closed outcome anticipated by the plan's risk section.

### Potential issues

No implementation defects remain. The campaign emits existing RFI/NaN warnings during
some prepared-spectrum variants; they are deterministic, do not abort a calculation,
and the exact isolated rerun proves they do not destabilize the serialized result.

## Manual Testing Required

No manual testing remains.

- `chime_rigorous_acf_fits.png`: all four central peaks, residual structures, total
  uncertainties, and excluded labels were inspected.
- `dsa_rigorous_acf_fits.png`: all four DSA peaks, structured residuals, fit curves,
  and excluded labels were inspected.
- `chromatica_modulation_qualification.png`: every conditional amplitude is shown as
  excluded on a log scale; the pathological totals remain visible.
- `chromatica_cross_band_fit.png`: all diagnostic centers are visible, no fit is drawn,
  and the no-joint-fit annotation matches the JSON.

Hashes and panel-level observations are recorded in the analysis review files.

## Recommendations

### Critical

None.

### Important

None for implementation readiness. Scientifically, do not quote the superseded
`alpha = 3.63 +/- 0.39` result or any current Chromatica modulation index.

### Nice to Have

- A future performance-only refactor could cache prepared window variants between
  reruns; it must preserve the current exact reproducibility and cache fingerprinting.

### Follow-Up Work

- Treat the failed gates as experiment-design evidence: improved dynamic-spectrum
  calibration or more independent bandwidth may be needed before another physical DSA
  modulation claim.
- Keep PR #200 draft and change its description from a provisional cross-band fit to a
  fail-closed campaign result unless a later, separately reviewed dataset qualifies.

## References

- [Plan](plan-chromatica-rigorous-scintillation-campaign.md)
- [Implementation](implement-chromatica-rigorous-scintillation-campaign.md)
- [Research](research-chromatica-rigorous-scintillation-campaign.md)
- [Experiment](experiment-chromatica-scintillation-fit-uncertainty.md)
