# Burst energetics — flux-calibration & cross-band validity review

Research note (codebase investigation, 2026-06-22). Scope: does the absence of absolute
flux calibration / per-band normalization invalidate the E_iso table produced by
`analysis/calculate_burst_energies.py`? Reviewer: Claude (session ac2a9913), at user request.

## Verdict (TL;DR)

| Question | Answer | Confidence |
|---|---|---|
| Q1. Are `c0_C` and `c0_D` in commensurable units? | **No** — independent per-telescope arbitrary units | High (code-confirmed) |
| Q2. Is summing `E_CHIME + E_DSA` defensible? | **No** — adds two independently-scaled bands; corrupts even cross-burst *ranking* of the total | High |
| Q3. Is the `γ_D ≈ −5` rail a DSA bandpass-rolloff artifact? | **Likely** — the `c0/γ` fit path applies no bandpass flattening | Medium (plausible mechanism; needs a raw-spectrum check) |

**Bottom line:** the table cannot be published as absolute energies, and the CHIME+DSA *sum*
is not even a reliable relative-ordering statistic. Salvageable only as a single
calibrated-band, order-of-magnitude estimate (see Recommendations).

## Methodology

Traced the amplitude `c0` from the joint-fit output back through the model likelihood to the
data ingest, and inspected the stored posteriors. Files: `burstfit_joint.py`, `burstfit.py`,
`analysis/scattering-refit-2026-06/run_joint_fit.py`, `joint_json/*_joint_fit.json`.

## Evidence

1. **`c0` is fit in per-channel S/N units (corrected 2026-06-22).** Earlier this note said `c0`
   "inherits whatever units the `.npy` carries." Direct read of the loader corrects that: before the
   data reaches `FRBModel`, `BurstDataset._bandpass_correct` (`io.py:131-137`) **z-scores every
   channel** by its own off-pulse mean/std, and `io.py:143-146` downsamples without renormalizing,
   keeping "units as S/N." So the fitted `c0` is the per-channel **signal-to-noise** spectral
   amplitude, with the arbitrary instrumental gain divided out — not raw native units. This does
   **not** weaken the verdict: S/N is still telescope-specific (the S/N→Jy factor σ_S = SEFD/(√(n_pol·Δν·Δt)·G)
   differs per band per channel), so `c0_C` and `c0_D` remain incommensurable and unsummable. It
   sharpens the fix: multiply each band's S/N by its own σ_S(ν) (the radiometer step) — see
   `docs/rse/specs/plan-radiometer-flux-cal.md`.

2. **No flux/units metadata anywhere.** `*_joint_fit.json` keys are
   `[burst, alpha, tau_1ghz, log_evidence, log_evidence_err, alpha_bounds, percentiles, ncall]`
   — no `units`, `calibration`, or `Jy` field. Nothing asserts the `.npy` is in Jy·ms.

3. **`c0_C` and `c0_D` show no common scale** (medians):

   ```
   burst        c0_C    c0_D        burst        c0_C    c0_D
   chromatica  0.719   0.107        oran        0.827   0.082
   freya       0.743   0.294        phineas     0.342   0.361
   hamilton    0.124   0.049        whitney     0.077   0.253
   isha        1.451   0.027        wilhelm     0.127   0.763
   johndoeII   0.131   0.219        zach        1.929   0.185
   mahi        0.059   4.181
   ```
   Each band is ~O(0.01–4) in its own normalization, set by its own telescope's pipeline. CHIME
   (transit cylinder) and DSA-110 (separate array) share no flux standard — there is no
   cross-band tie anywhere in the config or loader.

4. **The pipeline's own design treats the per-channel amplitude/spectrum as a nuisance.** The
   gain-marginalized path (`burstfit_joint.py:95-97`, `burstfit.py:728-742`) integrates the
   per-channel gain out analytically — "the gain absorbs the burst spectrum AND scintillation…
   c0 and gamma drop out of the sampled vector." The energetics calc resurrects `c0/γ` from the
   *other* (sampled) path and treats them as physical flux densities. This is an internal
   contradiction: the pipeline marginalizes the spectrum away precisely because the per-channel
   amplitude is not trusted as an absolute physical quantity.

5. **No bandpass flattening in the `c0/γ` path.** A smooth power law is fit to the channel
   amplitudes as-is. Any uncorrected frequency-dependent gain (band-edge rolloff) is absorbed
   into `γ`. The DSA `γ_D` posteriors rail at the −5 prior floor for several bursts (Phineas
   median −4.998, posterior `[-4.9995, -4.9954]`, err ±0.001-0.003; Chromatica/Zach/Hamilton/Oran
   similar), and Whitney sits at the opposite extreme (`γ_D=+2.91`). A railed/extreme index that
   tracks the band edge is the signature of an uncorrected bandpass, not an astrophysical spectrum.

## Findings per question

**Q1 — commensurability.** `c0_C` and `c0_D` are amplitudes in independent, uncalibrated,
per-telescope units (Evidence 1–3). They are not comparable; their ratio carries no physical
meaning.

**Q2 — the sum.** `E_iso = E_CHIME + E_DSA` adds `c0_C·(CHIME scale)` to `c0_D·(DSA scale)` with
two *different* unknown scales. This is worse than "absolute scale unknown": because the
CHIME:DSA energy mix swings per burst (Isha 0.5% DSA → Wilhelm 74% DSA), the unknown per-band
factors reshuffle the bursts differently, so even the *relative ordering* of the total is not
preserved. The code's caveat ("a common unknown factor leaves ordering intact") assumes ONE
factor; there are two.

**Q3 — the rail.** The `c0/γ` path does not flatten the bandpass, so an uncorrected DSA rolloff
would be absorbed into `γ_D` and drive it to the prior floor — a plausible, possibly dominant,
explanation for the −5 rail. Confirm by overplotting a raw DSA channel-amplitude spectrum against
the fitted power law (or by comparing to the gain-marginal `gain_spectrum`, which is the
calibration-robust per-channel estimate).

## Recommendations (prioritized)

1. **Do not sum the two bands.** Drop `E_iso = E_CHIME + E_DSA` from the table.
2. **Anchor one band and report it alone.** CHIME publishes catalog fluences (Jy·ms, with
   realistic beam-model uncertainties). Tie each CHIME waterfall to its catalog fluence, then
   report **CHIME-band `E_iso`** as an order-of-magnitude value with explicit uncertainty.
3. **Fix the bandpass before trusting `γ`.** Flatten the per-channel gain (or use the
   gain-marginal `gain_spectrum`) and re-fit; the −5 rail may relax. Until then, mark railed
   `γ_D` as prior-limited, not measured.
4. **If neither band is calibratable:** do not report absolute energies. Report only
   calibration-robust quantities (τ, Δν_d, spectral *shape* with caveats) or fluence *ratios*
   within a single band.
5. **Add posterior-propagated error bars** to whatever energy is reported; for railed parameters
   report a one-sided limit.

## Unblock progress (2026-06-22, beam/SEFD hunt)

**DSA beam — found, real cube.** `DSA110_beam_1.h5` (measured Jones E-field, X/Y feeds,
41 freq × 1801 θ × 73 φ; freq axis 1.2–1.6 GHz despite a wrong `freq_Hz`/`MHz` label). Loader +
trilinear `G(θ,φ,ν)` in `analysis/dsa_beam.py` (boresight=1; gain(1.8°,1.4 GHz)=0.477 = half power,
matching the 3.6° phi-avg FWHM). Replaces the analytic Airy fallback in
`dsa110-continuum/.../beam_model.py`. Path seam: `telescopes.yaml beam_model_h5`.

**DSA SEFD — measured per-epoch, not a constant.** Canonical source is the **dsa110-rt** SEFD
dashboard (`github.com/dsa110/dsa110-rt`, served on `lxd110h23:5777`): per-day × per-calibrator ×
baseline-bin (`sefd_0-200m`…`sefd_800-2500m`), median ~few-thousand Jy (warn >5000). Consistent
with the continuum repo's measured `5800 Jy/element` (T_sys=25 K, N_ant=96, η=0.7); the sim
config's "300–500 Jy" is a stale placeholder. **Action:** pull the measured median SEFD nearest
each burst's epoch from the h23 dashboard.

**CHIME — NOT calibrated upstream (verified 2026-06-22).** Direct inspection of a real
`singlebeam_*.h5` via the canonical `baseband_analysis` reader (h17, see
`docs/rse/specs/research-chime-singlebeam-flux-units.md`) **disproves** the earlier assumption
that the singlebeam product is already in Jy. The singlebeam is a `BBData` object
(`tiedbeam_baseband`/`tiedbeam_power`) with **no `units` attribute** and only a complex per-input
*timing/beamforming* gain applied (`calibrator = gain_…taua_ref_cyga_timing.h5`; the whole effect
is `core/calibration.py:66` `baseband *= conj(gain)`, coherent-summation phasing, no SEFD, no
primary-beam, no Jy). `tiedbeam_power` median ≈ 6.5×10⁷ — arbitrary correlator units. The
singlebeam pipeline (`pipelines/form_singlebeam.py`) has **no flux step**. So CHIME needs the
*same* radiometer machinery as DSA (Andersen+2023 / Michilli+2021 baseband-at-known-position:
primary-beam model at the burst position + system SEFD), run as a separate downstream step.
`flux_jy_per_unit_C` stays `null`.

## Remaining to verify

- Per-burst DSA **pointing offsets** (θ,φ from beam centre) to evaluate `G` per sightline.
- ~~CHIME `singlebeam_*.h5` flux units (confirm Jy)~~ — **done**: not Jy, pre-flux-cal
  (`docs/rse/specs/research-chime-singlebeam-flux-units.md`). Now need the CHIME primary-beam model
  + per-event SEFD to *apply* the flux cal.
- ~~Whether the FLITS input `.npy` preserved the power scale or renormalized to S/N~~ — **resolved
  2026-06-22**: `BurstDataset` z-scores each channel to per-channel **S/N** (`io.py:131-146`), for
  both bands. So the eventual calibration is S/N × σ_S(ν) (radiometer), per
  `docs/rse/specs/plan-radiometer-flux-cal.md`.
- Direct bandpass check: raw DSA spectrum vs fitted `γ_D` for Phineas (the worst 30%-DSA case).

## Prior-art recon — how the field calibrates these bursts (verified refs)

Sweep 2026-06-22 (perplexity scholar + ADS/SciX; every reference below resolved via the ADS API,
BibTeX in `references.bib`). The path to a quotable energy is standard radio practice, not a
backend dig.

**The one equation, both bands.** Signal-to-noise → flux density via the radiometer equation:

```
S_nu = (S/N) * SEFD / [ sqrt(n_pol * d_nu * d_t) * G(theta,phi,nu) ]
SEFD = 2 k_B T_sys / A_eff          fluence F = integral S_nu dt  [Jy*ms]
```

The only missing ingredients are each telescope's **SEFD** and the **beam response `G`** at the
burst's sky position. FLITS already whitens by off-pulse MAD, so the data is effectively in S/N
units — it is one SEFD multiply (× beam correction) away from Jy.

**CHIME (400–800 MHz) — templates:**
- Andersen+ 2023 (AJ 166, 138), *Flux Calibration of CHIME/FRB Intensity Data* — the dedicated
  method paper (radiometer + primary-beam model; the equation above).
- CHIME/FRB Collab. 2024 (ApJ 969, 145), *Updating the First CHIME/FRB Catalog with Baseband
  Data* — recomputes Catalog-1 fluences from **lower limits** to true values once the baseband
  position fixes the beam gain. Catalog-1 numbers are lower limits precisely because position
  uncertainty leaves `G` unknown.
- Michilli+ 2021 (ApJ 910, 147), baseband pipeline — flux cal at a **known** position against
  steady calibrators (Crab, Cas A). Co-detections have good positions, so this regime (not the
  real-time intensity one) applies.

**DSA-110 (1.3–1.5 GHz) — templates:**
- Law+ 2024 (ApJ 967, 29), *DSA Science: First FRB and Host Galaxy Catalog* — DSA's catalog/
  flux-cal paper; SEFD, calibrators (3C286/3C48) + noise-diode scale, single-burst fluence.
  **Gap:** the standalone DSA-110 *system* paper (Ravi et al.) is still unpublished, so the
  SEFD/beam description lives in Law+2024 (+ Connor+ 2023, ApJ 949, L26) and the DSA-110 thesis —
  cite those, not an instrument paper.

**Energetics — the field does NOT sum raw bands.** Standard band-limited isotropic energy
(Zhang 2018, ApJ 867, L21):

```
E_iso = 4*pi*D_L(z)^2 / (1+z) * F_nu * d_nu          (flat-spectrum form; the 1/(1+z) is the
                                                       bandwidth k-correction)
```

For a power law F ∝ nu^beta the k-correction generalizes to
`K = integral_rest nu^beta d_nu / integral_obs nu^beta d_nu` with nu_rest = (1+z) nu_obs.
Best practice for two bands is a **joint spectral fit on a common absolute (Jy) scale**, then one
integrated `E_iso` — *not* two independently-scaled per-band energies added together. FLITS
already does a joint fit; it just needs calibrated inputs. (z-from-DM via Macquart+ 2020,
Nature 581, 391, is not needed here — these sightlines have spectroscopic z.)

**Two concrete errors this surfaces in the current `burst_energies.tex`:** (i) it applies **no
`(1+z)`** k-correction; (ii) it **sums** two uncalibrated bands (compounding the Q2 problem).

## Related (separate lane, noted not reviewed)

A concurrent agent reworked the foreground/DM-budget search (`scratch/codetection/why_missed.py`;
`DEFAULT_CLUSTER_IMPACT_KPC = 5000` split from the 100 kpc galaxy gate). Physics is sound
(clusters contribute DM at Mpc impact; galaxies need <100 kpc for scattering) and the user
validated it — out of scope for this energetics review.
