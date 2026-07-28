# CHIME factor-ladder ACF analysis contract

**Status:** software-qualified interface; no burst result or science claim

**Factors:** 16, 32, 64, 128, 256, 512

**Execution gate:** the Zach pilot and 72-product campaign remain blocked on
ratified per-burst dispersion measures, owner-reviewed masks, product
generation, and owner review.

## One analysis interface

`scintillation.scint_analysis.acf_fitting` owns:

- factor and frequency-grid validation;
- equal-channel, equal-signal-to-noise, and fixed sub-band plans;
- single and multiple Lorentzian fits in physical frequency lag;
- additive versus multiplicative model labels;
- modulation-index algebra and a noise-debiased direct estimator.

`fit_lorentzian_components` accepts an ACF, frequency lags in MHz, optional
per-lag errors, channel width, and factor. It fits positive lags only. Lag zero
is excluded because its radiometer-noise spike is not scintillation. All width
bounds and starting points derive from the physical lag grid, not a factor's
bin number.

For a factor-tagged product, pass both:

```python
result = fit_lorentzian_components(
    lags_mhz,
    acf,
    acf_err=acf_err,
    max_components=3,
    channel_width_mhz=0.390625 / upchannel_factor,
    upchannel_factor=upchannel_factor,
)
```

The returned record includes factor, channel width, fit domain, model-selection
criterion, component widths and errors, and modulation parameterization. Its
Bayesian-information-criterion plus approximate F-test component count is
diagnostic only because per-lag error bars do not encode covariance between
neighboring ACF lags.

For a science-facing model comparison, `acf_evidence.compare_acf_evidence`
uses a correlated-lag covariance and nested-sampling evidence. It accepts the
same factor tag and validates the same grid. Its physical two-screen candidate
uses the multiplicative cross-term model. The nested sampler is an optional
dependency and must be present in the qualified execution environment.

## Sub-band definitions

`build_subband_plan` returns contiguous channel slices plus physical center
frequencies, bandwidths, and allocation weights.

- `equal_channels`: balanced channel counts. For a fixed factor, physical
  bandwidths are nearly equal.
- `equal_snr`: balances
  `sum((max(signal, 0) / noise_rms)^2)`. Signal is the on-pulse spectrum after
  off-pulse subtraction. Noise is the per-channel uncertainty of that
  on-minus-off mean, using the measured off-pulse standard deviation and exact
  valid-sample/time weights. This calculation assumes independent time samples;
  a measured temporal covariance requires a covariance-aware replacement.
  Masked, invalid, or zero-noise channels contribute zero.
- `fixed`: reuses complete, contiguous, owner-reviewed channel intervals.

The earlier splitter accumulated signed signal and called the result equal
signal-to-noise. That was only an approximation. The new plan records the exact
weight definition. Sub-band boundaries must be recomputed for every factor;
channel indices from one factor are not portable to another.

Time windows likewise remain physical-time inputs. Convert them to samples
separately for every factor because the upchannelizer trades frequency
resolution for time resolution.

## Lorentzian and modulation parameterization

Every component is

```text
m_i^2 / (1 + (lag / bandwidth_i)^2)
```

so the fitted non-negative `m_i` is the square root of that component's
extrapolated zero-lag ACF excess.

The selected multiple-Lorentzian fit is a **phenomenological sum**:

```text
baseline + sum_i m_i^2 rho_i
```

Its zero-lag equivalent is `sqrt(sum_i m_i^2)`. This is useful for describing
multiple ACF scales, but it is not automatically a physical multiple-screen
modulation index. Its uncertainty is propagated from the full fitted
component-modulation covariance, including correlations.

For independent multiplicative screens, use the separately named model:

```text
ACF + 1 = product_i (1 + m_i^2 rho_i)
```

It retains every cross term. For two fully modulated screens, the zero-lag ACF
excess is 3 and the total modulation index is `sqrt(3)`, not `sqrt(2)`.
The evidence result derives the total modulation for every posterior sample,
then reports its weighted median and 16th/84th-percentile uncertainties.
Choosing that physical model for burst inference remains a science decision;
the ordinary component fitter does not select it automatically. The separate
evidence path compares it against the single-screen model only after its
injection-calibrated science gate is available.

Factor changes do not change the Lorentzian equation. They do change
resolution-dependent limits. The deterministic fitter uses the observed lag
spacing and fit span for its width bounds. The evidence fit records its
factor-dependent prior explicitly: lower width `0.5 * channel_width`, upper
width `0.25 * fitted_bandwidth`. Cross-factor comparisons must therefore
compare posterior support and censoring under the recorded prior, not only
best-fit numbers.

The direct estimator must subtract noise variance:

```text
m = sqrt(max(observed variance - noise variance, 0)) / mean signal
```

A non-positive intrinsic variance is a non-detection, not `m = 0`. Existing
plain `std/mean` outputs remain explicitly labeled diagnostic because they
include radiometer noise and baseline dilution. When uncertainties on the
mean and both variances are supplied, the direct estimator propagates them.

## Qualification

Synthetic known-truth tests cover:

- all six factor grids and factor-independent recovery of one physical width;
- rejection of unsupported factors and mismatched channel widths;
- exact grid coverage for equal-channel and fixed plans;
- equal total signal-to-noise squared for the noise-aware plan;
- recovery and selection of one and two separated Lorentzian scales;
- the additive `sqrt(2)` versus multiplicative `sqrt(3)` two-screen case;
- noise-variance subtraction and the non-detection branch.

These tests qualify software behavior only. They do not validate any
per-burst factor, mask, time window, ACF scale, modulation index, or scientific
interpretation.
