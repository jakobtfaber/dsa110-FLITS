# Chromatica cross-band scintillation fit

This analysis combines the four accepted `chromatica_hi` CHIME/FRB narrow-component
bandwidths with a fresh DSA-110 ACF/Lorentzian pass. It does **not** read the legacy
DSA `stored_fits` block.

The DSA fit selects four equal-S/N subbands, but the 1459.62 MHz narrow component is
excluded from the cross-band regression because its off-pulse null fails. The three
other DSA narrow components pass the component, off-pulse, and low-lag gates.

## Reproduce

From the FLITS repository root, with the locked environment installed:

```bash
NUMBA_DISABLE_JIT=1 .venv/bin/python \
  analysis/scintillation-dsa-lorentzian-2026-07-07/run_dsa_lorentzian_fits.py \
  --flits-root "$PWD" \
  --output-dir analysis/chromatica-cross-band-scintillation-2026-07-17/dsa \
  --band dsa --bursts chromatica --max-components 3

.venv/bin/python \
  analysis/chromatica-cross-band-scintillation-2026-07-17/fit_cross_band.py \
  --chime-json analysis/window-tuning-campaign-2026-07-17/results/chromatica_hi_campaign.json \
  --dsa-json analysis/chromatica-cross-band-scintillation-2026-07-17/dsa/chromatica_dsa_lorentzian_fits.json \
  --output-dir analysis/chromatica-cross-band-scintillation-2026-07-17/results

.venv/bin/python -m pytest -q \
  analysis/chromatica-cross-band-scintillation-2026-07-17/test_fit_cross_band.py
```

The result JSON records SHA-256 hashes for both analysis inputs, the raw DSA `.npz`,
and the fitter. The primary
model includes intrinsic log scatter because the formal no-extra-scatter power law
fails its chi-square goodness-of-fit test. The formal fit is retained in full; it is
not silently promoted as an acceptable exact law.
