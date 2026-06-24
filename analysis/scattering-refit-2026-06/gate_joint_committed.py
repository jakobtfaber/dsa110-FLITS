"""Gate the committed joint CHIME-DSA scattering fits against the runtime 3-level
FLITS fit-quality contract, reusing the authoritative classify_fit_quality.

Reads joint_json/{burst}_joint_fit.json + the paired {burst}_joint_ppc.json,
applies Level-1 physical bounds + prior-rail detection, Level-2 reduced-chi2
(classify_fit_quality, worst of the two bands), and Level-3 alpha-physics, then
writes joint_gate_verdicts.{csv,md} and a per-burst {burst}_joint_gate.json
(mirroring the *_fit_results.json field shape). These are standalone verdict
artifacts -- the fit-verify workflow globs *_fit_results.json, so it does NOT
auto-discover these *_joint_gate.json files. tau x dnu (Level-3) is not evaluable
here -- no per-sightline scintillation bandwidth -- so it is reported N/A.

Level-1 here is the physical-bounds + prior-rail subset of the contract; the
optimizer-convergence / Jacobian-conditioning Level-1 gates do not map to these
nested-sampling fits (no Jacobian; convergence is proxied by log_evidence_err).
"""

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # repo root
from scattering.scat_analysis.burstfit import classify_fit_quality  # noqa: E402

ALPHA_MIN, ALPHA_MAX = (
    1.0,
    6.0,
)  # L1 hard gate: alpha<1.0 achromatic (tau prop nu^-a meaningless); ADR-0004
ALPHA_SUBKOLM = 2.0  # 1.0<=alpha<2.0 = sub-Kolmogorov: L3 MARGINAL (inspect), not L1 FAIL; ADR-0004
TAU_MIN, TAU_MAX = 1e-4, 100.0  # ms
RAIL_EDGE = 0.1  # alpha within this (absolute) of a prior bound => prior-railed (unconstrained)
RAIL_SIGMA = 3.0  # OR median within this many sigma of a prior bound => rail-MARGINAL regardless of value; ADR-0004
KOLM_LO, KOLM_HI = 3.5, 4.5  # Level-3 PASS-consistent alpha window (Kolmogorov ref 4.0)
_RANK = {"FAIL": 2, "MARGINAL": 1, "PASS": 0}


def _worst(*flags):
    return max(flags, key=lambda f: _RANK[f])


def gate_one(burst, fit, ppc):
    """Classify one joint fit. ppc may be None (chi2 unknown -> Level-2 MARGINAL)."""
    alpha = fit["alpha"]["median"]
    tau = fit["tau_1ghz"]["median"]
    lo, hi = fit["alpha_bounds"]
    sig_lo = fit["alpha"].get("err_minus") or 0.0
    sig_hi = fit["alpha"].get("err_plus") or 0.0
    # railed if pinned in absolute terms OR the median sits within RAIL_SIGMA*sigma of either
    # prior bound -- a wide posterior reaching a bound is unconstrained, not a measurement (ADR-0004)
    rail = (
        min(alpha - lo, hi - alpha) < RAIL_EDGE
        or (alpha - lo) <= RAIL_SIGMA * sig_lo
        or (hi - alpha) <= RAIL_SIGMA * sig_hi
    )

    # Level 1 -- physical bounds (any failure => FAIL, regardless of chi2)
    l1_fail = []
    # ADR-0004: hard FAIL only below ALPHA_MIN (achromatic) or at/above the ceiling;
    # 1.0 <= alpha < 2.0 is admitted (sub-Kolmogorov) and handled as L3 MARGINAL below.
    if alpha < ALPHA_MIN or alpha >= ALPHA_MAX:
        l1_fail.append(f"alpha={alpha:.3f} outside [{ALPHA_MIN},{ALPHA_MAX})")
    if not (TAU_MIN < tau < TAU_MAX):
        l1_fail.append(f"tau={tau:.4g}ms outside ({TAU_MIN},{TAU_MAX})")
    l1 = "FAIL" if l1_fail else "PASS"

    # Level 2 -- reduced chi2 per band via the runtime classifier, worst of two.
    # A missing PPC, or a present-but-incomplete one (a chi2_* key absent), is
    # treated the same: chi2 unknown -> MARGINAL, rather than crashing the run.
    cc = ppc.get("chi2_chime") if ppc else None
    cd = ppc.get("chi2_dsa") if ppc else None
    if cc is None or cd is None:
        note = "no PPC (chi2 unknown)" if ppc is None else "incomplete PPC (chi2 missing)"
        l2, l2_note = "MARGINAL", note
    else:
        fc = classify_fit_quality(cc)[0]
        fd = classify_fit_quality(cd)[0]
        l2 = _worst(fc, fd)
        l2_note = f"chi2_C={cc:.2f}({fc}) chi2_D={cd:.2f}({fd})"

    # Level 3 -- alpha physics consistency. Physical bounds are Level 1's job
    # (don't re-impose a second, stricter alpha cut here); Level 3 only judges
    # closeness to Kolmogorov: in-window => PASS-consistent, else MARGINAL
    # (off-Kolmogorov). tau x dnu not evaluable here (no per-sightline dnu_d).
    l3 = "PASS" if KOLM_LO <= alpha <= KOLM_HI else "MARGINAL"

    # A prior-railed alpha is unconstrained -- the fit pinned the prior bound
    # rather than measuring alpha -- so it is never better than MARGINAL even if
    # every level otherwise passes (consistent with excluding the rail-pinned
    # chromatica/freya/hamilton joint fits as non-measurements). And because the
    # Level-3 tau x dnu check is not evaluable here (no per-sightline dnu_d), the
    # contract is only partially verified, so the ceiling is capped at MARGINAL --
    # a fit is never certified PASS on an incomplete contract check.
    final = _worst(l1, l2, l3, "MARGINAL" if rail else "PASS", "MARGINAL")
    if l1_fail:
        reason = "L1 " + "; ".join(l1_fail)
    elif l2 == "FAIL":
        reason = "L2 catastrophic " + l2_note
    else:  # final is MARGINAL (capped); list whatever drove it beyond the cap
        bits = []
        if rail:
            bits.append(f"alpha prior-railed (within {RAIL_EDGE} of bound) -> unconstrained")
        if l2 == "MARGINAL":
            bits.append("L2 " + l2_note)
        if l3 == "MARGINAL":
            tag = "sub-Kolmogorov (inspect)" if alpha < ALPHA_SUBKOLM else "off Kolmogorov"
            bits.append(f"L3 alpha={alpha:.2f} {tag}")
        if not bits:  # passed every evaluable level; MARGINAL only because of the cap
            bits.append("L3 tau x dnu not evaluable (no dnu_d) -> capped at MARGINAL")
        reason = "; ".join(bits)

    return {
        "burst": burst,
        "alpha": alpha,
        "tau": tau,
        "rail": rail,
        "chi2_chime": cc,
        "chi2_dsa": cd,
        "l1": l1,
        "l2": l2,
        "l3": l3,
        "tau_dnu": "N/A (no dnu_d)",
        "final": final,
        "reason": reason,
    }


def _is_allexp(fit):
    """ADR-0003 PBF family: legacy mixed-PBF fits predate pbf_C/pbf_D and omit them."""
    return fit.get("pbf_C") == "exp" and fit.get("pbf_D") == "exp"


def main():
    import os

    jdir = Path(__file__).resolve().parent / "joint_json"
    fits = {fp: json.loads(fp.read_text()) for fp in sorted(jdir.glob("*_joint_fit.json"))}
    # ADR-0003 fail-closed: the committed joint_json fits are mixed-PBF, whose alpha are superseded.
    # Regenerating citable verdicts from them would leak a superseded sub-Kolmogorov alpha (e.g. oran,
    # whitney) into the energy/citable path -- calculate_burst_energies reads these *_joint_gate.json
    # quality flags and drops only FAIL. Refuse to overwrite the committed verdicts from mixed input
    # unless the caller explicitly opts into a clearly NON-CITABLE interim run.
    mixed = [fp.name for fp, fit in fits.items() if not _is_allexp(fit)]
    if mixed and os.environ.get("FLITS_GATE_ALLOW_MIXED") != "1":
        print(
            f"!! {len(mixed)}/{len(fits)} joint_json fit(s) are mixed-PBF (pbf_C/pbf_D != 'exp'); their"
            "\n   alpha are superseded (ADR-0003). Refusing to overwrite the committed citable verdicts"
            "\n   from mixed input. Re-grade from the canonical all-exp fits (pending the canonical-fit-"
            "\n   set decision), or set FLITS_GATE_ALLOW_MIXED=1 to write NON-CITABLE interim verdicts.",
            file=sys.stderr,
        )
        return 1
    provenance = "mixed-PBF-interim (NON-CITABLE)" if mixed else "all-exp"
    rows = []
    for fp, fit in fits.items():
        burst = fp.name.replace("_joint_fit.json", "")
        ppc_fp = jdir / f"{burst}_joint_ppc.json"
        ppc = json.loads(ppc_fp.read_text()) if ppc_fp.exists() else None
        v = gate_one(burst, fit, ppc)
        rows.append(v)
        worst_chi2 = None
        if v["chi2_chime"] is not None and v["chi2_dsa"] is not None:
            worst_chi2 = max(v["chi2_chime"], v["chi2_dsa"])
        (jdir / f"{burst}_joint_gate.json").write_text(
            json.dumps(
                {
                    "burst": burst,
                    "pbf_provenance": provenance,
                    "chi2_reduced": worst_chi2,
                    "chi2_chime": v["chi2_chime"],
                    "chi2_dsa": v["chi2_dsa"],
                    "tau": v["tau"],
                    "alpha": v["alpha"],
                    "alpha_railed": v["rail"],
                    "quality_flag": v["final"],
                    "notes": [v["reason"]],
                },
                indent=2,
            )
        )

    cols = [
        "burst",
        "alpha",
        "rail",
        "tau",
        "chi2_chime",
        "chi2_dsa",
        "l1",
        "l2",
        "l3",
        "tau_dnu",
        "final",
        "reason",
    ]
    out = Path(__file__).resolve().parent
    with (out / "joint_gate_verdicts.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    md = [
        "# Joint-fit quality verdicts (committed " + jdir.name + ")",
        "",
        f"**PBF provenance: {provenance}.** alpha floor = {ALPHA_MIN} (ADR-0004; 1.0<=a<2.0 = sub-Kolmogorov MARGINAL).",
        "",
        "| " + " | ".join(cols) + " |",
        "|" + "---|" * len(cols),
    ]
    for r in rows:
        md.append(
            "| "
            + " | ".join(f"{r[c]:.3f}" if isinstance(r[c], float) else str(r[c]) for c in cols)
            + " |"
        )
    n = {f: sum(r["final"] == f for r in rows) for f in ("PASS", "MARGINAL", "FAIL")}
    md += [
        "",
        f"**{len(rows)} fits**: {n['PASS']} PASS / {n['MARGINAL']} MARGINAL / {n['FAIL']} FAIL.",
    ]
    (out / "joint_gate_verdicts.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    sys.exit(main())
