# Chromatica rigorous cross-band scintillation campaign

This analysis now applies one post-ACF measurement contract to the CHIME/FRB and
DSA-110 Chromatica products. Telescope-specific preprocessing and CHIME harmonic
masking remain separate, but both bands use the weighted one-to-three Lorentzian
selector, moving-block bootstrap with model reselection, fit/window variants,
alternative-shape check, off-pulse and low-lag controls, matched-data injection, and
the same fail-closed reporting schema.

No legacy DSA `stored_fits` value enters the result. The matched injections add a
known thin-screen source spectrum to real off-pulse spectra of the same duration,
preserving the observed channel mask, channelization, radiometer noise, and production
ACF normalization.

## Result

The rigorous campaign admits zero of four DSA-110 subbands:

- 1321 MHz: window/fit-policy width instability and failed matched recovery.
- 1351 MHz: model-reselecting bootstrap width interval is too broad.
- 1396 MHz: window/fit-policy width instability and failed matched recovery.
- 1460 MHz: off-pulse null, alternative-shape, and matched-recovery failures.

It also admits zero of four CHIME/FRB subbands under the same stricter contract. Two
fail matched-injection calibration; two fail declared width/shape variants. Component
count stability is recorded separately and continues to block modulation indices even
where a narrow width might otherwise be stable.

Consequently, `results/chromatica_cross_band_fit.json` records the joint fit as
unavailable and reports no cross-band exponent. The earlier provisional
`alpha = 3.63 +/- 0.39` characterization is superseded and must not be quoted.
No `m_narrow`, `m_broad`, or `m_total` value is admitted for either telescope.

## Reproduce

From the FLITS repository root:

```bash
NUMBA_DISABLE_JIT=1 uv run --frozen python \
  analysis/chromatica-cross-band-scintillation-2026-07-17/run_rigorous_campaign.py \
  --root "$PWD" \
  --output-dir analysis/chromatica-cross-band-scintillation-2026-07-17/rigorous

uv run --frozen python \
  analysis/chromatica-cross-band-scintillation-2026-07-17/fit_cross_band.py \
  --campaign-json \
    analysis/chromatica-cross-band-scintillation-2026-07-17/rigorous/chromatica_rigorous_scintillation.json \
  --output-dir analysis/chromatica-cross-band-scintillation-2026-07-17/results

uv run --frozen pytest -q \
  scintillation/scint_analysis/tests/test_rigorous_campaign.py \
  scintillation/scint_analysis/tests/test_window_campaign.py \
  analysis/chromatica-cross-band-scintillation-2026-07-17/test_fit_cross_band.py
```

The campaign JSON records the random seeds, trial counts, thresholds, code/input
hashes, package versions, per-gate reasons, full central ACFs, compact variant fits,
and the exact commands. Figures were reviewed against
`docs/dev/figure-review-protocol.md`; their hashes and panel-level verdicts are in the
two `figures.review.json` files.
