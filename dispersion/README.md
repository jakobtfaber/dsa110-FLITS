# dispersion/ — status of the DM estimators

The manuscript's adopted DMs come from the phase-coherence campaign
preserved at `Faber2026/analysis/dm-joint-phase-v2/` (adoption record:
`Faber2026/analysis/docs/rse/specs/verified-dm-adoption-2026-07-13.md`).
Nothing in this package is the source of a quoted manuscript DM.

Modules kept here because current code imports them (superseded as
estimators, load-bearing as libraries):

- `chime_dm.py` — arrival-time regression estimator (superseded); its
  `K_DM` constant and `exgauss` kernel are imported by the v2 estimator
  and the injection tooling.
- `dm_power_analysis.py` — in-tree DM-power variant (superseded); its
  manifest loader (`load_manifest_rows`) is imported by the v2 campaign
  runner.
- `dmphasev2.py` — in-tree DM-phase variant (superseded); its
  `DMPhaseEstimator` is imported by
  `scattering/scat_analysis/dm_preprocessing.py`.
- `dm_phase_analysis.py` — battery-era harness (superseded); imported by
  kept regression tests (`tests/test_dmphase_recovery.py`,
  `tests/test_dm_power.py`).

The battery campaign package (`dm_campaign/`), the vendored published
packages (`external/DM_phase`, `external/DM-power`), and their tests are
archived at `.archive/outdated-science/dm-battery-2026-07/`.
