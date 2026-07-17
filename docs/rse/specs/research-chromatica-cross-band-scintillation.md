# Research: Chromatica cross-band scintillation fit

Date: 2026-07-17

## Question

Can the accepted CHIME `chromatica_hi` scintillation bandwidths and a fresh DSA-110
measurement be described by one power law,

\[
\gamma(\nu)=\gamma_{1\,\mathrm{GHz}}(\nu/1\,\mathrm{GHz})^\alpha?
\]

## Evidence inspected

- CHIME source: `analysis/window-tuning-campaign-2026-07-17/results/chromatica_hi_campaign.json`.
  It contains four accepted, resolved narrow-component measurements at 623.6--749.2 MHz.
- DSA source data: `scintillation/data/chromatica.npz` (external-data symlink).
- Fresh DSA driver:
  `analysis/scintillation-dsa-lorentzian-2026-07-07/run_dsa_lorentzian_fits.py`.
  The driver explicitly removes legacy YAML `stored_fits`, recomputes ACFs from the
  `.npz`, evaluates 2, 3, and 4 equal-S/N subband candidates, and selects the largest
  viable candidate.

## Findings

The fresh DSA run selects four subbands. The narrowest component in every subband
passes the component-level width/modulation flags and the low-lag-excision test.
However, the 1459.62 MHz width fails the off-pulse null: the on-pulse width is
1.38918 MHz and the off-pulse median is 1.38549 MHz. That scale is therefore
instrumental and must not enter the cross-band regression.

The remaining DSA measurements are:

| Frequency (MHz) | gamma (MHz) | fit error (MHz) | off-pulse null | low-lag stability |
|---:|---:|---:|:---:|:---:|
| 1321.063 | 0.72846 | 0.08979 | pass | pass |
| 1351.097 | 0.59571 | 0.08014 | pass | pass |
| 1395.889 | 2.06527 | 0.15769 | pass | pass |
| 1459.620 | 1.38918 | 0.10084 | **fail** | pass |

The DSA driver's burst-level `measurement_status` remains `measurement` for
non-CHIME data even when an off-pulse subband fails. The cross-band fitter therefore
must enforce the per-subband null and stability results itself.

For CHIME, the regression uncertainty should include the reported fit, finite-scintle,
and window-selection terms in quadrature. For DSA, the current fresh driver exposes
the Lorentzian covariance error and the artifact-control decisions; it does not expose
an analogous window-campaign systematic.

A no-extra-scatter weighted log-linear fit to the seven admitted points has a very
poor goodness of fit (chi-square about 68 for 5 degrees of freedom). A single exact
power law is therefore rejected under the stated measurement errors. The primary
reported cross-band characterization should include a fitted log-scatter term, while
also preserving the rejected no-scatter fit and leave-one-DSA-out sensitivities.

## Decision

Implement a fail-closed, provenance-recording fitter that:

1. admits only accepted/resolved CHIME narrow components;
2. combines all three CHIME uncertainty terms in quadrature;
3. admits only DSA `gamma_1` components with no quality flags, a passing off-pulse
   null, and passing low-lag stability;
4. reports both the formal weighted power law and a power law with intrinsic log
   scatter;
5. labels the intrinsic-scatter result as the primary characterization because the
   formal no-scatter model fails its goodness-of-fit test; and
6. records every admitted and rejected point, input hash, fit diagnostic, and
   sensitivity result.
