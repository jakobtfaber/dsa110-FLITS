# Audit: CHIME-side structure-DM extraction (retraction + rebuild)

**Date:** 2026-06-24
**Trigger:** expert visual review of the DM-phase diagnostics — bursts are not vertical, the
sweep is large, the "DM-phase" curves are flat/oscillatory, and it is universal across all 12.
**Status:** Pillar-2 DMs RETRACTED on main (nulled, `dm_confidence=under-audit`); rebuild pending.
Pillar-4 (positions) UNAFFECTED.

## What was wrong (confirmed by `scripts/diag_dedisp.py` on zach, in docker)

singlebeam: `dt=2.56 µs`, `ntime=55949` (~143 ms; pre-dedispersed at record time, burst at ~77 ms),
871 ch over 400.4–799.2 MHz.

1. **Inter-channel unit bug (1000×).** The extraction rolled channels by
   `1e-3 * K_DM * (1/f² − 1/f_ref²) * dm_c` with `K_DM=4148.81`, f in MHz. The physical full-band
   delay @DM=262 is **5.086 s = 1,986,647 samples**; the recipe rolled by **5.086 ms = 1987 samples**
   — **ratio exactly 1000.0**. The `1e-3` factor is wrong for f in MHz. Damage was partly masked
   because the data is already pre-dedispersed into a 143 ms window, but the band is still mis-aligned.
2. **DM-phase recovered nothing.** Every curve is near-flat/oscillatory (flat_ratio 1.2–1.5); the
   robust-argmax just picks the global max of noise wiggles, and the bootstrap σ is small only because
   the noise pattern is stable across resamples. A real structure-DM peak is sharp (ratio ≫ 1). This
   holds even on zach, where proper dedispersion shows a clear **SNR~15** burst.
3. **Downstream signal-loss.** The RFI-mask → coarse-peak (savgol argmax) → ±512 window → time-flip →
   z-score-display chain throws away a burst that `coherent_dedisp(bb, dm_c)` (proper, `time_shift=True`)
   renders clearly. Diagnostic: proper-dedisp collapse peak SNR~14.8 vs recipe-path SNR~12.1, but the
   extraction figures show noise.
4. **The recovery test was blind to all of this.** `tests/test_dmphase_recovery.py` dispersed *and*
   de-dispersed with the same `1e-3*K_DM` constant, so it recovered the injected DM without ever
   checking absolute physical scaling. It must be rebuilt to disperse with the *physical* constant.

## Rebuild plan (full rebuild on library dedispersion — user-approved scope)

Re-derive on `baseband_analysis` proper dedispersion; do not hand-roll inter-channel alignment.

- **P1 — physical-scale recovery test (test-first).** Rewrite `test_dmphase_recovery.py` to disperse a
  synthetic burst with the PHYSICAL constant (`t[s]=4.148808e3*DM*(f_MHz^-2 − f_ref^-2)`) and assert the
  estimator recovers it. Make the same convention the single source of truth for the estimator.
- **P2 — dedisperse correctly.** Use `coherent_dedisp(bb, dm_c)` (default `time_shift=True`), validated
  against zach (must reproduce the SNR~15 vertical burst). Search a small *physical* residual-DM grid for
  the structure-max; no `1e-3` factor, no manual roll.
- **P3 — make the curve peak.** On a known-signal burst (zach), confirm the DM-phase curve has a sharp
  interior peak (flat_ratio ≫ the ~1.3 noise floor) before trusting any other burst. Gate every burst on
  curve sharpness, not just "interior".
- **P4 — re-extract all 12 + figure-review from scratch + re-validate**; only then un-null the DMs that
  genuinely peak, and re-assess which co-detections yield a real CHIME structure-DM.

## Expert methodology verdict (2026-06-24, astronomy-astrophysics-expert) — METHOD PIVOT
After P1 (bug fix) + P2 prototype: with proper `coherent_dedisp` and the fixed estimator, zach's
structure-max curve is flat_ratio ≈ 1.97–2.27 across windows ±4/±12/±30 ms (synthetic sharp-structure
control gives 48.5). Verdict:
- **flat_ratio ~2 = NON-detection of structure, NOT a wrong DM.** Structure-max (DM_phase; Seymour/
  Michilli/Hessels+2019) measures DM curvature only from unresolved temporal sub-structure; a smooth,
  scattered, low-S/N (~3–11) blob has none, so the curve is flat by construction even at the correct DM.
  Structure-max is the WRONG primary tool for these CHIME singlebeam bursts.
- **Circularity is real:** dedispersing at DM_DSA then fitting a tiny reference-centered residual (≈0) is
  mostly a windowing artifact, not an independent confirmation. The defensible statement is an EXCLUSION
  from a WIDE, reference-independent DM search, not "δ_DM≈0".
- **Right estimator is already in the repo:** `burstfit.py` M2/M3 forward-fit (DM `delta_dm` + scattering
  `tau_1ghz` fitted jointly), which is scattering-aware (removes the S/N-DM scattering bias) and works on
  smooth bursts. Use a wide flat DM prior (±20–50 pc/cm³), report the marginalised DM posterior, σ =
  σ_stat ⊕ σ_scat (DM shift τ-free vs τ-fixed). Run the repo PASS/MARGINAL/FAIL gates.
- **Per-burst outcome, not global:** bursts whose posterior is tight enough → report DM±σ + |ΔDM| 95%
  exclusion (support Pillar 2). Bursts too low-S/N → "CHIME singlebeam does not independently constrain
  DM" (lean on Pillar 4 position). Keep structure-max only as a secondary cross-check gated on a
  null-based ≥5σ peak significance (phase-scramble / off-pulse null), NOT a flat_ratio number.
- **Caveat to confirm:** what DM the baseband was coherently dedispersed at vs what we re-dedisperse at —
  `delta_dm` is only meaningful relative to a known reference; intra-channel coherent smearing at the
  wrong DM is not undone by an incoherent re-trial.

### Revised plan (supersedes P2–P4 above)
- **P2′** — drop structure-max as primary. Prototype a `burstfit` M2/M3 fit of the CHIME band-collapsed
  (or 2-D) burst for zach over a WIDE reference-independent DM prior; confirm a marginalised DM posterior
  + validation PASS.
- **P3′** — define the agreement test (pre-registered): |DM_CHIME−DM_DSA|/σ_CHIME < 2 AND a stated |ΔDM|
  95% exclusion; classify each burst constrains / does-not-constrain.
- **P4′** — fit all 12, validate each, report the honest split; un-null only the bursts that constrain DM.

## Custom DM tool (user directive 2026-06-24: "rewrite into a custom tool that works 100%")
Built `dispersion/chime_dm.py` — a clean-room, telescope-agnostic, pure numpy/scipy DM estimator
(no flits deps → runs in the baseband docker image and host). Method = the textbook scattering-aware
DM measurement, NOT structure-max:
1. **wide incoherent DM search** (reference-independent) → band-collapsed peak-S/N coarse DM that can
   be far from `dm_ref` (so a real offset is FOUND → de-circularises the agreement test);
2. **scatter-corrected arrival-time regression**: per sub-band fit an exponentially-modified Gaussian
   (Gaussian ⊗ one-sided exp(−t/τ)); the Gaussian centre `t0` is the scattering-DECONVOLVED arrival
   time → weighted linear fit of `t0` vs `K_DM·(ν⁻²−ν_ref⁻²)`; slope = residual DM, covariance×χ²_red
   = honest σ_DM. Few sub-bands / large σ → `constrains_dm=False` (no fabricated value).

**Validation (`tests/test_chime_dm.py`, 12 tests, independent numerical injector):** known-DM recovery
in BOTH CHIME (400–800) and DSA (1281–1531) bands across τ∈{0,2,5 ms}; wide-search recovers a +30
offset (de-circularisation); non-detection floor flagged at S/N~1; DSA σ_DM > CHIME σ_DM (narrow-band
lever arm); σ calibrated (pull scatter < 3). **Real DSA data:** phineas_dsa (6144×5121, DM 610.274) →
S/N 140, 7/8 sub-bands, **DM = 610.206 ± 0.006, constrains_dm=True** — recovers the catalogue DM.

**σ caveat (expert):** the reported σ is statistical (χ²_red-inflated). The 0.07 pc/cm³ phineas offset
is physically negligible but ~10σ on σ_stat alone — so the **agreement test (in association.py) must
apply a physical tolerance floor (~1 pc/cm³)**, not the raw pull. Tool reports the honest statistical
measurement; the floor lives in the pillar-2 agreement logic.

### Remaining (P4″)
Docker-extract all 12 CHIME singlebeam bursts with `chime_dm.measure_dm` (coherent_dedisp at the DSA DM
→ measure), classify constrains/does-not-constrain, wire the agreement test with the physical floor,
regenerate the report, un-null only the bursts that constrain DM, PR to main.

## Provenance
Prior (retracted) numbers live in git (PR #29, commit cc64b7b) and off-repo at
`/data/.../results/chime_side_inputs.json`. Diagnostic: `scripts/diag_dedisp.py`,
`diagnostics/diag_dedisp_zach.png`.
