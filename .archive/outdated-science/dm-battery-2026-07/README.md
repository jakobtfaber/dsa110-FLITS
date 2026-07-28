# Archived: DM multi-method battery campaign (2026-07, superseded)

Historical code from the DM re-measurement battery. No manuscript DM traces
to this code. The adopted implementation is the phase-coherence campaign
preserved at `Faber2026/analysis/dm-joint-phase-v2/` (source branch
`agent/dm-phase-v2`, commit `c07f1f1`); the adoption record is
`Faber2026/analysis/docs/rse/specs/verified-dm-adoption-2026-07-13.md`,
and the archived campaign documents live in
`Faber2026/analysis/docs/rse/archive/dm/`.

Contents:

- `dm_campaign/` — injection harness, estimator adapters, battery and
  adaptive-arrival runners (plan: `plan-dm-measurement-methods.md`, archived
  as above). Its role was served: the injection matrix proved the earlier
  in-tree DM-power null on DSA was an implementation artifact.
- `external/` — vendored published packages `DM_phase` (GPL-3,
  `b7cf5fd61436`) and `DM-power` (`f7787355ca28`, no published license —
  internal reproduction only), consumed only by the battery adapters. See
  `external/README.md` for the vendor patches.
- `tests/` — the battery-era test files (`test_dm_injection.py`,
  `test_dm_adapters.py`, `test_adaptive_arrival.py`). They import
  `dispersion.dm_campaign` and will not run from here; retained for the
  record, not for collection.

Still live in `dispersion/` (deliberately NOT archived — load-bearing for
current code): `chime_dm.py` (K_DM used by the v2 estimator),
`dm_power_analysis.py` (manifest loader used by the v2 runner),
`dmphasev2.py` (used by `scattering/scat_analysis/dm_preprocessing.py`),
and `dm_phase_analysis.py` (imported by kept regression tests).

Historical result products (gitignored, never in this repo's history) were
moved to `~/Data/Faber2026/results-library/dispersion/archive/`.
