# Research: Chromatica rigorous scintillation campaign

**Date:** 2026-07-17
**Scope:** Internal codebase
**Code state:** `a0b3029`
**Related Documents:**
- [Chromatica cross-band research](research-chromatica-cross-band-scintillation.md)
- [Chromatica cross-band validation](validation-chromatica-cross-band-scintillation.md)

## Question / Scope

Determine whether the accepted CHIME/FRB and provisional DSA-110 scintillation
bandwidths and modulation indices were produced by the same measurement pipeline,
identify the gaps that prevent a physical DSA modulation-index claim, and define the
requirements for a common, fail-closed remeasurement campaign.

## Codebase Findings

### Shared ACF core, different end-to-end pipelines

Both paths ultimately use `calculate_acfs_for_subbands` and `calculate_acf`. The ACF
is normalized by `(mean_on - mean_off)^2` when an off-pulse level is supplied, and
its diagonal error combines a product standard error with a heuristic finite-scintle
term (`scintillation/scint_analysis/analysis.py:246-340`,
`scintillation/scint_analysis/analysis.py:555-680`). Equal-S/N subband construction is
also shared (`scintillation/scint_analysis/analysis.py:619-648`).

The accepted CHIME window campaign then uses a custom unweighted `curve_fit` fitter
over positive lags up to 5 MHz, with CHIME harmonic-comb masking, one- versus
two-Lorentzian BIC selection, a four-to-one scale-separation requirement, a
line-versus-Lorentzian shape gate, amplitude rails, and a resolution gate
(`scintillation/scint_analysis/window_refit.py:112-235`). It separately adds a
finite-scintle width uncertainty and excludes unphysical narrow-component amplitudes
from the frequency-law fit (`scintillation/scint_analysis/window_refit.py:319-390`).

The fresh DSA driver uses `ScintillationAnalysis` for preparation, disables the
synthetic-noise template and joint 2D fit, and applies the weighted
`compare_lorentzian_components` selector over a configurable fit span
(`analysis/scintillation-dsa-lorentzian-2026-07-07/run_dsa_lorentzian_fits.py:79-103`,
`analysis/scintillation-dsa-lorentzian-2026-07-07/run_dsa_lorentzian_fits.py:1338-1425`).
That selector fits one through three Lorentzians plus a constant using diagonal ACF
errors and requires both delta-BIC and an approximate nested F-test before adding a
component (`scintillation/scint_analysis/revalidation.py:277-425`).

Thus the ACF estimator is shared, but preprocessing, fit weighting, model families,
selection rules, fit spans, and uncertainty budgets are not.

### The DSA status gate is not fail-closed

The DSA driver calculates per-subband off-pulse and low-lag controls, but passes them
to a finalizer designed specifically for CHIME
(`analysis/scintillation-dsa-lorentzian-2026-07-07/run_dsa_lorentzian_fits.py:1407-1425`,
`analysis/scintillation-dsa-lorentzian-2026-07-07/run_dsa_lorentzian_fits.py:1518-1555`).
`finalize_measurement_status` explicitly returns the provenance status unchanged for
non-CHIME data, ignoring failed physical controls
(`scintillation/scint_analysis/chime_artifact_guards.py:480-513`). That is why the
Chromatica DSA burst remained labelled `measurement` while its 1459.62 MHz subband
failed the off-pulse null.

### Current modulation indices are conditional fit amplitudes

Both fitted models parameterize a component as
`m^2 / (1 + (lag/gamma)^2)`. A fitted amplitude can represent a physical modulation
index only if the ACF has a valid off-level normalization; the code documents that
without off-level subtraction the amplitude is noise-floor diluted
(`scintillation/scint_analysis/freya_scintillation.py:62-76`).

The accepted CHIME subbands are all two-component fits. The campaign's stored `m` is
the narrow-component amplitude, not the quadrature total. The DSA 1351 MHz fit also
selects two components: `m_narrow = 0.7605` and `m_broad = 1.0182`, but the broad
width exceeds its fit window and is flagged. Neither the narrow amplitude nor the
formal quadrature total can be quoted as an unconditional physical modulation index
when the decomposition is unstable.

The reported `m_err` values are conditional covariance errors from the selected
model (`scintillation/scint_analysis/revalidation.py:364-390`). They omit variation
from on/off windows, subband boundaries, retained-lag policy, alternative model
families, correlated ACF lags, and model selection. The DSA pass also lacks a
window/subband systematic analogous to the CHIME campaign's `gamma_win_sys`.

## Synthesis

A common campaign must share the fit and reporting contract while retaining only
physically necessary telescope-specific preparation. It must evaluate window,
subband, fit-span, and first-lag variants; apply per-subband off-pulse and low-lag
gates to every telescope; compare plausible ACF shape families; and propagate
variant plus resampling uncertainty into both bandwidth and amplitude.

The primary bandwidth should remain the narrowest admitted Lorentzian scale only when
that component is stable across the model and analysis variants. Modulation reporting
must distinguish `m_narrow`, `m_broad`, and `m_total`; each must have its own admission
status. A failed or unidentifiable broad component must not silently contaminate a
total modulation index.

The current PR #200 cross-band fit should remain diagnostic until the remeasurement
campaign either qualifies stable DSA measurements or fails closed. Failure to
qualify a DSA modulation index is a scientifically valid campaign result.

## References / Sources

- `scintillation/scint_analysis/analysis.py:246-340,555-680`
- `scintillation/scint_analysis/window_refit.py:112-235,319-390`
- `scintillation/scint_analysis/revalidation.py:277-425`
- `scintillation/scint_analysis/chime_artifact_guards.py:480-513`
- `scintillation/scint_analysis/freya_scintillation.py:62-76`
- `analysis/scintillation-dsa-lorentzian-2026-07-07/run_dsa_lorentzian_fits.py:79-103,1338-1555`
