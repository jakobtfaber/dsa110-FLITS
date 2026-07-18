# Experiment: Chromatica scintillation fit and uncertainty architecture

**Date:** 2026-07-17
**Author:** AI Assistant
**Status:** Complete
**Related Documents:**
- [Research: rigorous scintillation campaign](research-chromatica-rigorous-scintillation-campaign.md)
- [Plan: rigorous scintillation campaign](plan-chromatica-rigorous-scintillation-campaign.md)

## Experiment Goal

Choose a common central estimator and uncertainty method for a rigorous CHIME/FRB
and DSA-110 Chromatica scintillation campaign.

**Primary Question:** Do the current DSA widths depend materially on the CHIME-style
unweighted fitter, and are the stored conditional covariance errors adequate once
correlated ACF residuals and model reselection are included?

**Context:** The existing CHIME and DSA paths share the ACF calculator but use
different fitting and selection rules. A common campaign needs one defensible central
estimator without mistaking diagonal fit covariance for the full uncertainty.

## Hypothesis

The weighted DSA selector and the CHIME-style unweighted selector should give similar
central values for well-behaved single-scale ACFs. Moving-block residual resampling
should produce wider intervals and expose unstable component counts.

**Success Criteria:**
- Central widths agree to within 15% for subbands that both methods fit.
- The resampling method reports model-selection frequencies as well as parameter
  quantiles.
- The chosen production method preserves the ACF's lag ordering and reselects the
  model on each resample.

## Approaches Tested

### Approach 1: weighted BIC plus F-test selector

The current DSA `compare_lorentzian_components` path fits one through three
Lorentzians plus a constant, weighted by the production ACF diagonal errors.

### Approach 2: CHIME-style unweighted one/two-Lorentzian selector

The `window_refit._fit_subband` path was applied to the same four DSA ACFs without a
CHIME harmonic mask. This isolates fit weighting and model-selection behavior from
telescope preprocessing.

### Approach 3: weighted selector plus moving-block residual bootstrap

The selected weighted model supplied the center curve. Contiguous eight-lag residual
blocks were resampled 80 times per subband with seed `20260717`; the one-to-three
component model was reselected on every draw.

## Experiment Setup

**Environment:** locked FLITS environment at git commit `a0b3029`.

**Test Data:** the four freshly calculated DSA-110 Chromatica equal-S/N subband ACFs
from `scintillation/data/chromatica.npz`, using the checked-in DSA preprocessing and
25 MHz fit span.

**Execution:**

```bash
NUMBA_DISABLE_JIT=1 .venv/bin/python /tmp/dsa_fit_architecture_experiment.py \
  > /tmp/dsa_fit_architecture_experiment.json
```

The prototype was deliberately kept outside the repository; the raw JSON is copied
below so the decision does not depend on a transient file.

## Results

| Center (MHz) | Weighted model | Weighted gamma (MHz) | Unweighted gamma (MHz) | Bootstrap gamma 16/50/84% (MHz) | Bootstrap model counts |
|---:|:---:|---:|---:|---:|:---|
| 1321.063 | 1L | 0.7285 | 0.7229 | 0.4662 / 0.6914 / 0.9107 | 1L: 76, 2L: 4 |
| 1351.097 | 2L | 0.5957 | 0.5875 | 0.4609 / 0.5959 / 0.7779 | 1L: 2, 2L: 78 |
| 1395.889 | 1L | 2.0653 | 1.8205 | 1.5045 / 1.9462 / 2.4343 | 1L: 76, 2L: 4 |
| 1459.620 | 2L | 1.3892 | 1.3908 | 1.2296 / 1.3850 / 1.6336 | 1L: 16, 2L: 64 |

The corresponding narrow-component modulation-index bootstrap intervals were
`0.694/0.788/0.847`, `0.713/0.763/0.831`, `0.940/0.998/1.081`, and
`1.432/1.481/1.533` in increasing frequency order.

The weighted and unweighted centers agree within about 1% in three subbands. The
1396 MHz width differs by about 12%, still within the predeclared 15% comparison
criterion but large enough to require an explicit model-family systematic. The
moving-block intervals are approximately two to three times wider than the stored
conditional width errors. The 1460 MHz component count is unstable in 20% of draws;
that point also fails the independent off-pulse null in the existing analysis.

## Comparison Matrix

| Criterion | Weighted selector | Unweighted selector | Weighted plus block bootstrap |
|---|---|---|---|
| Uses production ACF errors | Yes | No | Yes |
| Central-value stability | Good | Similar, one 12% shift | Same as weighted center |
| Correlated-lag uncertainty | No | No | Approximate, preserves local residual correlation |
| Model-selection uncertainty | No | No | Yes, via reselection counts |
| Production suitability | Central estimator only | Sensitivity check | Central estimator plus uncertainty layer |

## Key Insights

1. The optimizer mismatch is not the principal source of overconfidence.
2. Conditional covariance errors materially understate the uncertainty visible in
   correlated residuals.
3. Model count must be treated as an uncertain outcome, especially for the excluded
   1460 MHz subband and multi-component modulation totals.
4. The 1396 MHz point needs explicit fit-span and alternative-shape sensitivity; it
   is the main leverage point in the existing cross-band result.

## Recommendation

Use the weighted BIC plus corroborating F-test selector as the common central
estimator. Add a moving-block residual bootstrap with model reselection, fixed seeds,
and explicit model-count frequencies. Treat the CHIME-style unweighted fit and a
generalized-Lorentzian fit as sensitivity checks, not as the production center.

The bootstrap is not a replacement for window, subband, fit-span, or matched-data
injection tests. Those enter the final qualification and systematic budget
separately. A modulation total is eligible only when the same component structure is
stable and every contributing component is resolved and physical.

## Appendix: Raw Experiment Data

```json
{
  "seed": 20260717,
  "n_bootstrap": 80,
  "block_length": 8,
  "runtime_seconds": 96.4,
  "subbands": [
    {"center_mhz": 1321.063, "weighted": {"n": 1, "gamma": 0.728457, "gamma_err": 0.089793, "m": 0.759436}, "unweighted": {"n": 1, "gamma": 0.722916, "m": 0.773021}, "bootstrap": {"n_counts": {"1": 76, "2": 4}, "gamma_q16_q50_q84": [0.466210, 0.691400, 0.910730], "m_q16_q50_q84": [0.694080, 0.788010, 0.847460]}},
    {"center_mhz": 1351.097, "weighted": {"n": 2, "gamma": 0.595714, "gamma_err": 0.080145, "m": 0.760500}, "unweighted": {"n": 2, "gamma": 0.587453, "m": 0.762791}, "bootstrap": {"n_counts": {"1": 2, "2": 78}, "gamma_q16_q50_q84": [0.460855, 0.595893, 0.777948], "m_q16_q50_q84": [0.713023, 0.762899, 0.831408]}},
    {"center_mhz": 1395.889, "weighted": {"n": 1, "gamma": 2.065269, "gamma_err": 0.157686, "m": 0.989769}, "unweighted": {"n": 1, "gamma": 1.820485, "m": 0.999157}, "bootstrap": {"n_counts": {"1": 76, "2": 4}, "gamma_q16_q50_q84": [1.504475, 1.946249, 2.434337], "m_q16_q50_q84": [0.939612, 0.998346, 1.080938]}},
    {"center_mhz": 1459.620, "weighted": {"n": 2, "gamma": 1.389182, "gamma_err": 0.100841, "m": 1.482294}, "unweighted": {"n": 2, "gamma": 1.390842, "m": 1.474920}, "bootstrap": {"n_counts": {"1": 16, "2": 64}, "gamma_q16_q50_q84": [1.229612, 1.384543, 1.633589], "m_q16_q50_q84": [1.432163, 1.480979, 1.533480]}}
  ]
}
```
