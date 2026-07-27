# Component-count ladder audit — how counts were chosen, and the D4 gap

**Date:** 2026-07-17. **Scope:** joint CHIME+DSA scattering mass-refit, 12 bursts.
**Purpose:** the honest, citable answer to "how were the per-burst temporal
component counts (C_n D_m) chosen?" — intended for the paper's methods appendix,
not just internal notes.

## The finding: D4 was never on the table
Across every fit launched in this campaign, the DSA (D-band) component count
histogram is:

| D-count | runs |
|--------:|-----:|
| D1 | 17 |
| D2 | 13 |
| D3 |  6 |
| **D4** | **0** |

`--components-D 4` was **never launched for any burst**. The ladder ceiling was
D3 (zach, phineas) and C4 (hamilton, CHIME side). Consequently no burst could
ever have been *found* to need a fourth DSA component — the hypothesis was never
tested.

## How counts were actually assigned (per burst)
Only two bursts had a genuine evidence ladder that compared neighboring counts:

- **isha** — D1 vs D2 at fixed s2 → D1 selected (D2 lost).
- **whitney β is gain-marginalization-dominated (task #12 classification).** At the
  evidence-confirmed C2D2 count, β sweeps the ENTIRE prior range purely by changing
  the gain-prior variance s2, with τ swinging 27× — the data do not constrain the
  screen index against the gain model:

  | gain treatment | β | τ_1GHz (ms) | lnZ |
  |---|---:|---:|---:|
  | flat / profiled (production convention) | 3.025 (+0.035/−0.018, **floor**) | 0.077 | 20160.07 |
  | s2 = 100 (regularized gain) | 3.429 (**interior**) | 0.037 | 20860.02 |
  | s2 = 10 (tight gain prior) | 3.988 (**ceiling**) | 0.989 | 15320.13 |

  This is a **gain-systematic-dominated** classification, NOT a “β=3 screen”
  measurement. **Valid comparisons:** the two fixed-s2 fits share a proper prior
  normalization, so s2=100 vs s2=10 IS a valid Bayes factor — ΔlnZ = **+5539.9**
  decisively favoring the interior (β=3.43) over the ceiling. The flat/profiled-gain
  lnZ (20160.07) uses a different (improper) gain-prior normalization and is NOT
  comparable to the fixed-s2 values, so no flat-vs-s2 Bayes factor is formed. Bottom
  line: among comparable gain priors the evidence prefers the interior screen; the
  production flat-gain “floor” is a convention-dependent corner of a gain-degenerate
  posterior, not a physical index. hamilton (C4D1, flat-gain β=3.003, τ→0) is being
  probed the same way (s2=10/100 + C5D1/C4D2 neighbors); if its floor also melts under
  a regularized gain, the floor-rail class dissolves into a shared whitney+hamilton
  gain systematic.
- **zach**: the D4 collapse at s2=100 and in the profiled C2D4 (ΔlnZ = −2.3 vs C2D3,
  chi²_D 1.13→1.12, screen params byte-identical, the +2.06 ms member left with a
  +6.6σ residual spike in BOTH D3 and D4 — figure `figs_multicomp/zach_d3_d4_resid.png`)
  is consistent with the resolution limit at 131 µs, not with the member being absent.
  Per the binning lever, the count is settled by a fine-binning refit, not this ladder.
  **Binning-lever diagnosis (2026-07-18):** the 131 us (t4) DSA binning is NOT set by
  the trailing-window cap (WIN_TRAIL_CAP_MS 30->12 is a no-op on the chosen t-factor).
  It is set by common_window=True, which unions DSA with CHIME; CHIME at 400-800 MHz
  scatters ~150x more, so its window runs to the ~44 ms record end and drags DSA to a
  44 ms window, so t_floor forces DSA to t4. Fix #1: common_window=False so DSA is
  independent. BUT (caught by owner review) DSAs own peak-anchored robust window is then
  only ~2.2 ms and TRUNCATES the cluster: the initial pulse sits at the peak and the
  +2.06/+2.52/+3.01 ms members fall OUTSIDE, because the ~1 ms quiet gap between the
  initial and the cluster exceeds WIN_MAX_GAP_MS (1.0). A truncated window breaks the
  count test by construction. Fix #2: a band-aware ENVELOPE window keyed on the full
  >5-sigma component span + margin, applied to the DSA band ONLY; CHIME keeps its
  original tail-following window (binning unchanged at t64/163.8 us, tail not clipped).
  Verified deployed-driver prep: CHIME 163.8 us (unchanged), DSA 32.8 us (t1, 4x finer)
  with a 5.90 ms window spanning [-1.4,+4.5] ms -> all four candidate components in-window.
  STANDING GUARDRAIL (owner): any windowed count test must first prove the window contains
  every candidate component; a peak-anchored window that keys on the brightest pulse only
  is the same failure class as the count shortfall itself. Fine pair (C2D3_fine/C2D4_fine,
  s2=100, _fine-suffixed) running (jobs 133/134); D4 verdict = fine-pair delta-lnZ only (guardrail 1).
- **mode-trap lesson for the methods appendix**: at s2=100, a lower-count fit can get
  trapped in a secondary (steep-β / runaway-τ) mode, inflating the apparent evidence gain
  of the next count. Read a count "win" as real only when the screen parameters are
  continuous across the count step; otherwise it is a mode difference, not a component.

## Fine pair (jobs 133/134) INVALID — off-window ghost components (2026-07-18)

The C2D3_fine / C2D4_fine pair is **INVALID as a D3-vs-D4 count test**; its
ΔlnZ(D4−D3) = +3549.6 must NOT be published or consumed anywhere. Mode-check passed
(both β healthy, 3.980 / 3.974) but the visual vet (owner standing rule: no verdict
from numbers alone — vet `/tmp/zachfine_vet.png`, deck slide 10) kills a naive read:

- BOTH fits park D1 **outside** the [0, 5.9] ms fitted window (t0 ≈ −8.3 and −9.4 ms,
  window-relative). The component C2D4 ADDS is D2 at t0 ≈ −1.36 ± 0.5 ms — **also
  off-window**. These are ghost components: only their scattered tails (if anything)
  reach the data, so the gain marginal exploits them as per-channel baseline / pedestal
  degrees of freedom. The +3550 is evidence for a model with unphysical off-window
  structure, NOT for the owner's 4th cluster member.
- The owner's cluster (+2.06/+2.52/+3.01 ms peak-rel = 3.56/4.02/4.51 window-rel) is fit
  as ONE broad component in BOTH fits, with identical ±4–6σ residual wiggles across
  3.4–4.7 ms in both — i.e. the cluster is present and **unfit by both** D3 and D4. So
  the count test never actually engaged the physics it was meant to.

The earlier "the 4th component absorbed tail that D3 mis-fit as scattering, dropping τ
0.29→0.183" reading is WRONG: the added component PRECEDES the window, so only its own
scattered tail leaks in as a pedestal (which is what let τ drop). A decisive number
pointing in the owner-expected direction (D=4), produced by an unphysical model, is the
trap — not a resolution.

**Root cause — t0 prior untethered from the window.** `burstfit.build_priors` sets
`t0 = (init.t0 − 2·max(init.tau_1ghz, 10.0), init.t0 + 2·max(init.tau_1ghz, 10.0))`.
Since tau_1ghz (ms) is always < 10, the `max(…, 10.0)` floors the half-width to 10 ms,
giving a t0 prior of **init.t0 ± 20 ms** regardless of the fitted window (~5.9 ms here).
All joint prior specs (`_joint_prior_spec`, `_joint_prior_spec_gain`,
`_joint_prior_spec_gain_gp`, `_joint_prior_spec_gain_shared_zeta` in burstfit_joint.py)
take `pC["t0"]`/`pD["t0"]` straight from build_priors, so every joint fit inherits the
±20 ms t0 range and can place components ≳10 ms outside the window.

**STANDING GUARDRAIL (owner, 2026-07-18): t0 priors bounded to the fitted window.**
Every temporal component's t0 prior must lie within its band's fitted window
[window_lo, window_hi] (optionally a small margin, NOT ±20 ms). This is the DUAL of the
existing "window must contain every candidate component" rule: the window must contain
the candidates, AND the prior must not let components escape the window. Off-window
kernels hand the gain marginal spurious per-channel baseline freedom and manufacture
ΔlnZ that is not about component count.

**Audit trail:** the invalid pair JSONs are snapshotted at
`flits-runs/data/joint/_invalid_zachfine_offwindow_20260718/` before any overwrite.

**Remediation (in progress, slotted behind two-screen Stage 0):**
(a) bound the joint t0 priors to each band's fitted window (fix in burstfit_joint.py
    prior specs / build_priors; validated + landed via PR, not hot-edited under running
    jobs); (b) re-run C2D3 vs C2D4 with bounded t0 and the 4th component seeded in-window
    near the cluster (window-rel ≈ 3.6/4.0/4.5 ms); (c) whether owner D=4 (initial + THREE
    cluster members) makes the honest grid C2D4-vs-C2D5 is a scope question pending
    team-lead sign-off before launch. Task #10 remains OPEN.


## v2 re-run harvest (jobs 169–182) — 2026-07-19

**Prior-spec:** v2 windowed t0 (PR #205). Do not cross-compare lnZ with any pre-2026-07-19 v1 evidence.
**Products:** `~/flits-runs/data/joint/{oran,johndoeII,zach}_joint_fit_*` (mtime 14:55–17:53 PDT);
vet figures `~/flits-runs/data/joint/_v2_harvest_20260719/{v2_harvest_vet,zach_v2_ladder_vet}.png`.
**All 14 jobs RC=0** (169–176 oran/johndoeII; 177–182 zach fine C2D3/D4/D5 × s2{10,100}).

### oran C1D1 vs C2D1 (production ghost remediation)

| arm | C1D1 lnZ / β | C2D1 lnZ / β | ΔlnZ(C2−C1) | C2 t0 / ζ |
|---|---:|---:|---:|---|
| s2=10 | 13634.95 / 3.965 | 13625.93 / 3.806 | **−9.0** | 16.65 / 0.021 |
| s2=100 | 13603.01 / 3.988 | 13603.16 / 3.929 | +0.1 | 19.00 / **737** |

v1 ghost had t0_C1 ≈ −5.2 ms off-window. Under v2 every t0 is in-window; the extra CHIME component is either fluence-null (s2=10) or ζ-runaway (s2=100). **Verdict: DROP to C1D1.** Production TOA CHIME structure for oran is single-component under v2.

### johndoeII C1D2 vs C2D2

| arm | C1D2 lnZ / β | C2D2 lnZ / β | ΔlnZ(C2−C1) | C2 t0 / ζ |
|---|---:|---:|---:|---|
| s2=10 | 11725.69 / 3.799 | 11723.33 / 3.748 | **−2.4** | 7.01 / **90** |
| s2=100 | 11679.91 / 3.866 | 11678.68 / 3.519 | **−1.2** | 11.12 / **757** |

v1 ghost had t0_C1 ≈ −6.2 ms, fluence 26 vs real 8. Under v2 the extra CHIME component is ζ-runaway on both arms; C1 wins both. **Verdict: DROP to C1D2.**

### zach fine ladder C2D3 / C2D4 / C2D5 (task #10)

| count | s2=10 lnZ / β | s2=100 lnZ / β |
|---|---:|---:|
| C2D3 | 25134.98 / **3.178** | 25074.60 / **3.165** |
| C2D4 | 26560.30 / **3.979** | 25064.50 / **3.166** |
| C2D5 | 26584.89 / **3.979** | 25081.73 / **3.161** |

**s2=10:** D3→D4 is a **MODE JUMP** (β 3.18→3.98, ΔlnZ +1425). INVALID as a count Bayes factor (same trap class as pre-clamp +3550, different mechanism: ceiling-β family vs interior). D5−D4 within the ceiling family is only +25.

**s2=100 (mode-continuous β≈3.16):**
- ΔlnZ(D4−D3) = **−10.1** → D3 preferred
- ΔlnZ(D5−D4) = +17.2; ΔlnZ(D5−D3) = +7.1
- Every count still carries one high-ζ null-like DSA member (ζ ∼ 170–450, broad t0 or dead fluence)
- Real DSA peaks under s2=100: initial ~1.5 ms + cluster ~3.5–4.0 ms; owner’s three-member cluster is not cleanly resolved as three tight components

**Verdict (task #10 under v2):** owner D=4 is **NOT supported**. Prefer **C2D3** on the mode-continuous s2=100 arm; do not publish s2=10 D4/D5 ΔlnZ. Residual structure around the cluster remains (see `zach_v2_ladder_vet.png`) but is not resolved by adding D4/D5 under windowed priors.

### Standing consequences

1. Production counts: oran → C1D1; johndoeII → C1D2 (both were C2 ghosts under v1).
2. TOA table (task #6): rewrite oran/johndoeII CHIME rows as single-component; flag prior production rows SUSPECT→SUPERSEDED by v2.
3. zach task #10: close as **D3 stands under v2 fine binning**; D4/D5 not evidence-backed once mode-check applied.
4. Prior-spec split remains: never mix v1/v2 lnZ.
