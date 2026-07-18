# Implementation Summary: Chromatica rigorous scintillation campaign

**Date:** 2026-07-17
**Author:** AI Assistant
**Status:** Complete
**Plan Reference:** [plan-chromatica-rigorous-scintillation-campaign.md](plan-chromatica-rigorous-scintillation-campaign.md)

## Overview

Implemented a common, fail-closed CHIME/FRB and DSA-110 scintillation measurement
campaign for Chromatica. The result supersedes the provisional cross-band fit: no
subband in either telescope passes every declared width gate, no modulation index is
admitted, and no cross-band exponent is reported.

## Plan Adherence

The implementation follows the approved architecture: weighted one-to-three
Lorentzian central fits, moving-block residual bootstrap with model reselection,
window/fit-policy variants, alternative narrow-component shape checks, per-subband
off-pulse and low-lag controls, and fixed-seed matched-data injection.

Two refinements were made during execution:

- Width stability and component-count stability are recorded separately. A stable
  narrow width need not be discarded solely because a broad component appears or
  disappears, but every modulation value remains ineligible unless the decomposition
  is stable.
- Injection coverage is evaluated against the final declared width interval, not the
  old conditional covariance. The latter was the quantity the experiment established
  as underestimated.

Both refinements make the implemented gates match the physical quantities they
qualify; neither changes the predeclared recovery-bias or coverage thresholds.

## Phases Completed

### Phase 1: Common fit and report contract

- Added `rigorous_campaign.py` with ordered component semantics, separate width and
  modulation eligibility, generalized narrow-shape sensitivity with fixed broad
  terms, block bootstrap, named uncertainty terms, and explicit-True qualification.
- Added eight focused unit tests, including known-truth Lorentzian recovery and a
  two-scale alternative-shape regression.

### Phase 2: Telescope adapters and variants

- Extended `window_refit.refit` additively to expose JSON-ready ACF arrays and optional
  prepared-spectrum objects, while retaining all existing return fields.
- Added the Chromatica campaign runner with fixed-boundary pointwise variants and
  campaign-level 2/3/4 equal-S/N partition records.

### Phase 3: Matched-data injection and products

- Exposed the existing thin-screen spectrum simulator through a documented public
  wrapper.
- Added known sources to real off-pulse spectra of the same subband and duration,
  retaining real channel masks, channelization, baseline/noise, and ACF normalization.
- Generated central ACF, modulation eligibility, manifest, and review products.

### Phase 4: Fail-closed cross-band regeneration

- Added a strict rigorous-schema loader; legacy per-telescope JSON cannot enter the
  production cross-band path.
- The cross-band builder now requires at least two qualified points per telescope and
  emits an explicit unavailable result otherwise.
- Regenerated JSON, CSV, PNG, and PDF outputs with no fitted exponent or envelope.

### Phase 5: Validation

- Focused and regression tests passed, including the ACF-evidence test under its
  declared `nested` extra.
- Ruff passed on new/touched conforming Python, and `git diff --check` passed.
- A fresh `uv run --isolated --frozen` campaign produced an exactly equal JSON result
  and byte-identical PNG figures.
- All four new PNG products were inspected at original resolution and reviewed.

## Files Created

- `scintillation/scint_analysis/rigorous_campaign.py`
- `scintillation/scint_analysis/tests/test_rigorous_campaign.py`
- `analysis/chromatica-cross-band-scintillation-2026-07-17/run_rigorous_campaign.py`
- `analysis/chromatica-cross-band-scintillation-2026-07-17/rigorous/*`
- `docs/rse/specs/research-chromatica-rigorous-scintillation-campaign.md`
- `docs/rse/specs/experiment-chromatica-scintillation-fit-uncertainty.md`
- `docs/rse/specs/plan-chromatica-rigorous-scintillation-campaign.md`

## Files Modified

- `scintillation/scint_analysis/acf_covariance.py`
- `scintillation/scint_analysis/window_refit.py`
- `scintillation/scint_analysis/tests/test_window_campaign.py`
- `analysis/chromatica-cross-band-scintillation-2026-07-17/fit_cross_band.py`
- `analysis/chromatica-cross-band-scintillation-2026-07-17/test_fit_cross_band.py`
- Cross-band result, README, and figure-review artifacts in the same analysis directory.

No files were deleted, and the parent `Faber2026` submodule pointer was not changed.

## Scientific Result

The DSA-110 qualification failures are:

| Center (MHz) | Central gamma (MHz) | Total sigma (MHz) | Failed gates |
|---:|---:|---:|---|
| 1321.063 | 0.6227 | 0.3777 | variant stability; matched injection |
| 1351.097 | 0.5198 | 0.4104 | bootstrap stability |
| 1395.889 | 1.9815 | 0.8071 | variant stability; matched injection |
| 1459.620 | 1.3892 | 0.3494 | off-pulse null; alternative shape; matched injection |

The values remain in the diagnostic record but are not measurements. In particular,
the 1351 MHz nominal center is not promoted merely because its median injection bias
is small: its model-reselecting width interval is too broad. The 1460 MHz result still
fails the independent off-pulse null.

All conditional `m_narrow`, `m_broad`, and `m_total` values are marked excluded. The
pathological large totals remain visible in the log-scale diagnostic so they cannot be
mistaken for omitted data.

## Reproducibility

The campaign JSON records seed `20260717`, 80 bootstrap draws per subband, 32 matched
injections per subband, block length 8, input and runner SHA-256 hashes, Python/NumPy/
SciPy versions, exact commands, thresholds, and all gate evidence. The isolated run
reported:

```text
campaign_json_exact_match=true
figure_png_exact_match=true
```

## Remaining Work

No implementation or scientific-validation work remains. Publication of the validated
commit and update of draft PR #200 are the final operational closeout steps.

## References

- [Research](research-chromatica-rigorous-scintillation-campaign.md)
- [Experiment](experiment-chromatica-scintillation-fit-uncertainty.md)
- [Plan](plan-chromatica-rigorous-scintillation-campaign.md)

**Implementation completed by AI Assistant on 2026-07-17.**
