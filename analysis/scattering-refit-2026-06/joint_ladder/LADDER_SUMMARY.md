# Joint CHIME–DSA scattering ladder — drain summary (2026-06-23)

All SLURM jobs drained (6 high rungs at 32 cores: 64528816–829; 8-core batch
903–915; all `COMPLETED`). Results pulled from
`hpcc:/central/scratch/jfaber/flits-runs/data/joint/` → `joint_ladder/`.

## Reading the evidence correctly

- **Profiled-s2 lnZ is empirical-Bayes (a profile Z). It is NOT comparable
  across component count N**, and the `(base)`/`sharedzeta` rows sit on a
  *different likelihood normalization* (~−10⁴) than the `CxDy` component
  runs (~+2×10⁴) — do not subtract across those families.
- **Component model selection uses the fixed-s2 ladder only.** A component is
  real iff ΔlnZ(N+1 vs N) is consistently positive (≳5) across s2 ∈ {1,10,100}.
  A sign flip with the prior scale = the component is prior-driven, not real.
- **Caveat (Codex review, 2026-06-23):** even `sharedzeta` vs `base` is *not* a
  pristine Bayes factor — they differ in likelihood class **and** α prior (e.g.
  casey base α∈[1,6] vs sharedzeta α∈[1.5,6]). Preferring sharedzeta is
  scientifically reasonable but the lnZ gap is not a clean factor.
- **PBF systematic (whole ladder, found via whitney 2026-06-24).** All joint fits
  used the script default **CHIME=powerlaw / DSA=exp** — physically incoherent
  (same burst, same sightline ⇒ one PBF shape; only τ scales with ν). For whitney
  the *physical* configs split hard: all-exp α=5.12 vs all-powerlaw α=1.51 (railed),
  and the data prefer **all-exponential by ΔlnZ ≈ +2708** (clean same-data/same-dim
  Bayes factor). PBF is therefore a **dominant** α systematic, not a detail — and
  the mixed default's extra (unphysical) shape freedom inflates lnZ (whitney mixed
  beats all-exp by +18), so its lnZ is *not* admissible. **Every ladder α below was
  obtained under the mixed default and needs an all-exp rerun before it is final.**
- **Only whitney / johndoeII / phineas(CHIME) have fixed-s2 grids.** Every other
  component verdict below (oran, isha, mahi, zach, phineas-DSA) rests on
  **profiled** cross-N Δ only — the very comparison flagged as not-clean above.
  Both bursts we *did* test with fixed s2 had their profiled hint **overturned**,
  so treat profiled-only component counts as PROVISIONAL pending an s2 grid.

## Fixed-s2 cross-N verdicts (the clean Bayes factors)

| burst | test | ΔlnZ across s2 = 1/10/100 | verdict |
|---|---|---|---|
| whitney   | C2D1 vs C1D1 | −5.3 / +1.7 / −4.6 | 2nd comp NOT real → single *(default prep; SUPERSEDED — see whitney note: under re-prep the 2nd comp is DSA, not CHIME, and is strongly real)* |
| whitney *(re-prep, all-exp)* | C2D2 vs C2D1 | +2706 / +2683 / +2671 | **2nd DSA comp REAL** — flat across s2, no sign flip (C2D1 rails 1.52 at every s2; C2D2 off-rail 5.65–5.68) |
| johndoeII | C2D2 vs C2D1 | +150.9 / −126.3 / +119.7 | **2nd DSA comp NOT real** (overturns profiled Δ+6.7) |
| phineas   | C3D1 vs C2D1 | +1.7 / +6.2 / −1.2 | 3rd CHIME comp not robust at D1 (DSA D2/D3 untested at fixed s2) |

## Per-burst result

α gate (validation contract Level 1): physical bound 1.5 < α < 6.0; **α at a
prior edge ⇒ FAIL** (railed posterior is not a measurement).

| burst | model | α | τ₁GHz (ms) | flag | note |
|---|---|---|---|---|---|
| **freya**     | sharedzeta | 4.36 ± 0.035 | 0.119 | PASS | near-Kolmogorov; PPC χ² 1.30/1.03 — model tracks both bands |
| **mahi**      | C1D1 (single; C2D1 Δ−291) | 3.80 +0.31/−0.24 | 0.278 | PASS | single comp preferred; overlay clean |
| **phineas**   | C3D3 (DSA multi-comp) | 3.33 ± 0.06 | 0.302 | PASS | α robust across ladder; overlay 1.02/1.34 |
| **chromatica**| sharedzeta | 3.29 ± 0.04 | 0.196 | PASS | PPC 1.14/1.16 (mild DSA tail under-fit) |
| **oran**      | C2D1 (Δ+178 **profiled only**) | 2.69 ± 0.16 | 0.733 | PASS (α) / PROVISIONAL (2-comp) | 2nd comp profiled Δ only — **no fixed-s2 grid**; overlay tracks |
| **wilhelm**   | sharedzeta | 2.56 ± 0.04 | 0.251 | **MARGINAL** | **PPC DSA χ²=4.55 (>3 fail)** — high-SNR shape mismatch; α a forced compromise |
| **zach**      | C2D3 (Δ+1451 vs C2D2) | 2.41 ± 0.02 | 0.682 | PASS | **profile-bias: 3.32→2.41**; but α s2-sensitive (2.76 at s2=1/10); D2 prior-dependent (fig-3) |
| **casey**     | sharedzeta | 2.40 ± 0.014 | 0.061 | PASS | PPC 1.41/0.99 |
| **isha**      | C2D1 (Δ+536 vs C1D1) | 5.39 +0.50/−1.95 | 0.017 | MARGINAL | α poorly constrained (DSA broad/low-SNR); τ very short |
| **hamilton**  | none usable | 1.50 (railed) | 0.045 | **FAIL (complete)** | every parameterization unphysical: base joint α=1.01 (below floor), sharedzeta & C1D1 rail at 1.50 |
| **whitney**   | C2D2, **all-exp PBF** (re-prep: 64-ch CHIME + 5σ DSA crop) | 5.12 +0.17/−0.16 | 0.109 | **PASS** (prep+PBF-corrected) | railing was *prep + PBF*. Default prep over-binned CHIME and buried the 2nd DSA peak → C2D1 railed α=1.52. Re-prep (finer CHIME + tight crop) reveals **2 DSA peaks** (t0 0.55/0.89 ms); C2D2 unrails α. PBF is decisive: **all-exp (thin screen) α=5.12** vs all-powerlaw α=1.51 (railed) — data prefer exp by **ΔlnZ +2708**. The mixed default (pw/exp) is unphysical and gave 5.21 (≈exp by luck). χ² 1.12/1.53, residuals clean. **Component reality confirmed:** fixed-s2 C2D2-vs-C2D1 ΔlnZ = +2706/+2683/+2671 (s2 1/10/100) — flat, no flip ⇒ 2nd DSA comp real. Caveat: differs from default-pipeline prep+PBF. |
| **johndoeII** | C2D1 (2nd DSA not robust) | 1.06–1.95 (unstable) | ~0.77 | **FAIL/MARGINAL** | every parameterization low (base 1.06, railed 1.50, best off-rail C2D2 1.95) — data genuinely want low α, not a pure model artifact |

## Headline

- **7 clean α (PASS, model-vs-data reviewed):** freya 4.36, mahi 3.80,
  phineas 3.33, chromatica 3.29, oran 2.69, zach 2.41, casey 2.40. (oran's
  *2-component* count is provisional — profiled-only, no fixed-s2 grid.)
- **+1 recovered via re-prep + PBF fix: whitney 5.12 +0.17/−0.16** (all-exp PBF;
  was complete FAIL). Recovery needed finer-CHIME + tight-DSA-crop **and** dropping
  the unphysical mixed PBF for all-exp — both must be applied campaign-wide before
  whitney (or any α) joins the others on equal footing. See whitney note + the PBF
  caveat above.
- **2 MARGINAL:** isha 5.39 (wide lower tail; DSA broad/low-SNR); **wilhelm 2.56
  (PPC DSA χ²=4.55 — α a forced compromise).**
- **Railing is parameterization-specific, not multi-vs-single** (C1D1 single
  rails too):
  - **hamilton — complete FAIL (figure-confirmed).** Overlay shows the **DSA band
    is essentially noise** (no real burst) and CHIME is one narrow spike → no
    scattering lever arm, so shared α rails 1.50. Every parameterization
    unphysical (base 1.01). No usable α.
  - **whitney — RECOVERED (prep + PBF corrected); was complete FAIL.** Under the
    default prep every parameterization railed (base α=3.75 flat, misses the CHIME
    multi-peak χ²=2.88; component path rails α=1.52). Two compounding causes, both
    *method* not burst: **(1) prep** — CHIME over-binned and the DSA burst sat in a
    wide noise window that buried its second peak; re-prepping (64-ch CHIME from
    native 1024; 5σ DSA crop) exposes **two DSA peaks**, and modeling them (C2D2)
    unrails α. **(2) PBF** — the default mixed CHIME=powerlaw/DSA=exp is unphysical;
    the physical configs give all-exp α=**5.12 +0.17/−0.16** (τ=0.109) vs
    all-powerlaw α=1.51 *railed*, and the data prefer **exp by ΔlnZ +2708**. Final:
    **α=5.12 (all-exp), C2D2**, χ² 1.12/1.53, residuals noise-like (figure-confirmed),
    logz +2700 over C2D1. A second profile-bias case, *opposite direction* from
    zach: the unmodeled component railed α **low** (→1.5); modeling it + the correct
    PBF pushed α **high** (5.12). Component reality **confirmed** by the fixed-s2
    grid: ΔlnZ(C2D2−C2D1) = +2706/+2683/+2671 at s2 1/10/100 (flat, no sign flip).
    Caveat: different prep+PBF from the rest of the ladder (the all-exp rerun fixes
    the PBF half campaign-wide).
  - **johndoeII — genuinely low/unstable.** All parameterizations low
    (1.06–1.95); the data want low α — not a pure model artifact. No clean number.
  - Where the component/shared-ζ path rails (hamilton, whitney, johndoeII), that
    path is mis-specified (candidate causes: cross-band DM/alignment error, or a
    band with no real scattering → degenerate lever arm) — investigate before
    quoting any component-model α.
- **zach** — *the 2.41 result is PBF-confounded (see all-exp section, 2026-06-24).*
  Under the mixed PBF, modeling hidden sub-components appeared to drop α from 3.32
  (single) to 2.41 (C2D3). But under the **physical all-exp PBF the C2D3 fit gives
  α=4.59 ± 0.04** (Δ+2.18) — so 2.41 was a mixed-PBF artifact, not a profile-bias
  measurement. **Settled (job 64539614):** the single-component all-exp α = **3.319
  ± 0.013** ≈ the single-component mixed 3.32, so the PBF is immaterial for the
  single-component fit but the multiplicity correction *reverses* sign — under the
  physical PBF modeling components *raises* α (3.319→4.59), the opposite of the
  mixed-PBF artifact (3.32→2.41). The "hidden component biases α high" thesis does
  not survive. **Rescue resolved — FAILS (fixed-s2 grid, jobs 64542330–45):** the
  3rd DSA component is prior-driven, not real — ΔlnZ(C2D3−C2D2) sign-flips across s²
  (+1443/−759/−0.4 at s²=1/10/100) and the C2D3 α is s²-unstable (2.40/2.77/2.76
  fixed vs 4.59 profiled). **zach = single-component all-exp α = 3.319 ± 0.013**;
  do not cite 2.41 or 4.59.
- **Component disagreements all resolved toward the Bayes factor:** oran/isha
  gain a component (evidence > visual); whitney drops to single and johndoeII's
  2nd DSA component is rejected (fixed-s2 sign-flip > profiled preference).

## All-exp PBF rerun — campaign PBF systematic (2026-06-24)

Whole ladder re-fit with the **physical all-exponential PBF** (each burst's
best-model, only the PBF changed from the mixed CHIME=powerlaw/DSA=exp default),
to isolate how much the unphysical mixed default biased α. Run as a workflow on
the cluster; α_exp vs α_mixed below.

| burst | α_mixed | α_exp (all-exp) | Δα | railed/wide? | note |
|---|---|---|---|---|---|
| whitney | 5.21 | 5.12 | −0.09 | no | re-prep+all-exp; trustworthy (done locally) |
| isha | 5.39 | 5.48 +0.42/−1.98 | +0.09 | **wide→upper** | unconstrained either PBF (MARGINAL stands) |
| freya | 4.36 | 4.356 ± 0.04 | −0.004 | no | **PBF-robust** |
| mahi | 3.80 | 3.17 +1.47/−1.18 | **−0.63** | **wide→floor** | α PBF-sensitive + poorly constrained under all-exp; re-check |
| phineas | 3.33 | 3.426 ± 0.05 | +0.096 | no | **PBF-robust** (job 64538743, 2h/16c rerun) |
| chromatica | 3.28 | 3.286 ± ~0.04 | +0.006 | no | **PBF-robust** |
| oran | 2.69 | 2.662 ± 0.16 | −0.028 | no | **PBF-robust** |
| wilhelm | 2.56 | 2.558 ± ~0.04 | −0.002 | no | **PBF-robust** (per-band PBF immaterial — α unchanged) |
| zach | 2.41 | **3.32 (single)** | — | no | C2D3 **rejected** (fixed-s2 sign-flip, jobs 64542330–45); single-comp all-exp α=3.319±0.013. 2.41 and profiled-C2D3 4.59 not citable |
| casey | 2.40 | 2.396 ± ~0.04 | −0.004 | no | **PBF-robust** |
| johndoeII | 1.58 | 1.573 (≈floor) | −0.007 | **railed** | hard-pinned at 1.5 under both PBFs |
| hamilton | 1.50 | 1.504 (≈floor) | +0.004 | **railed** | hard-pinned/unconstrained under both PBFs |

**Verdict (all 12 complete).** The PBF systematic is **small for the 7
well-constrained bursts** (freya, casey, chromatica, wilhelm, oran, whitney,
phineas): |Δα| ≤ 0.1, so their mixed-PBF α are confirmed, *not* overturned — the
all-exp values agree and become canonical for physical consistency. The
per-band-PBF preference (wilhelm ΔlnZ +4.0) is therefore **immaterial to the
science** (wilhelm α 2.56→2.558) and should be dropped as unphysical. **Five
block the ladder:** **zach** — the marquee profile-bias case — shifts **+2.18**
(2.41→4.59) under the physical PBF, so the published 2.41 ("hidden component biases
α low") was a **mixed-PBF artifact**; **mahi** (Δα −0.63, wide → floor); **isha**
(wide → upper); **johndoeII** & **hamilton** (hard-railed at 1.5, both PBFs). The
lesson: the PBF is negligible for clean sightlines but can dominate α for the
marginal ones (zach, mahi), so it must be physical before any α is quoted. zach's
single-component all-exp fit (job 64539614) gives α = 3.319 ≈ the single-component
mixed 3.32 — the PBF is immaterial for the single-component fit. The all-exp
fixed-s2 grid (jobs 64542330–45) then **resolved the rescue: it fails.** The 3rd
DSA component is prior-driven — ΔlnZ(C2D3−C2D2) sign-flips across s²
(+1443/−759/−0.4 at s²=1/10/100) and the C2D3 α is s²-unstable (2.40/2.77/2.76
fixed vs 4.59 profiled). **zach is therefore a single-component measurement,
α = 3.319 ± 0.013**; neither 2.41 (mixed) nor 4.59 (profiled C2D3) is citable.
