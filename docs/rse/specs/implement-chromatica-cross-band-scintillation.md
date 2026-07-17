# Implementation: Chromatica cross-band scintillation fit

Date: 2026-07-17
Plan: `docs/rse/specs/plan-chromatica-cross-band-scintillation.md`

## Implemented

- Locked a fresh Chromatica DSA-110 run under
  `analysis/chromatica-cross-band-scintillation-2026-07-17/dsa/`.
- Added `fit_cross_band.py`, which parses the CHIME campaign and fresh DSA result,
  enforces per-point admission gates, combines CHIME uncertainties, and fits both
  the formal and intrinsic-scatter power laws.
- Added per-point standardized residuals and leave-one-DSA-out sensitivity fits.
- Made near-zero intrinsic-scatter uncertainties fail closed as unidentifiable rather
  than serializing an overflow or infinity.
- Added a tested CLI and reproducibility README.
- Generated machine-readable point/result tables and title-free PDF/PNG figures.

## Scientific status encoded in the output

The DSA 1459.62 MHz point is retained in the point inventory but excluded because
its off-pulse null fails. The exact no-extra-scatter power law is retained but rejected
by its goodness-of-fit test. The intrinsic-scatter likelihood is designated as the
primary cross-band characterization.

## Verification commands

```bash
.venv/bin/python -m pytest -q \
  analysis/chromatica-cross-band-scintillation-2026-07-17/test_fit_cross_band.py

.venv/bin/python \
  analysis/chromatica-cross-band-scintillation-2026-07-17/fit_cross_band.py \
  --chime-json analysis/window-tuning-campaign-2026-07-17/results/chromatica_hi_campaign.json \
  --dsa-json analysis/chromatica-cross-band-scintillation-2026-07-17/dsa/chromatica_dsa_lorentzian_fits.json \
  --output-dir analysis/chromatica-cross-band-scintillation-2026-07-17/results
```
