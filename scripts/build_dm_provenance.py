#!/usr/bin/env python
"""Build the per-burst, per-telescope DM provenance table (V6 / Phase 6 P6.2).

Joins the CHIME-side DM measurements in ``crossmatching/chime_side_inputs.json``
with independent DSA-side DMs, documents the dedispersion method and producing
artifact for each side, and quantifies the CHIME-DSA agreement per burst
(``delta_dm`` and ``delta_dm_sigma``). Emits ``crossmatching/dm_provenance.csv``.

SOURCE CHANGE (2026-07-27): ``configs/bursts.yaml`` no longer carries the
independent DSA measurement — its dm/dm_err now mirror the adopted
CHIME-primary catalog (analysis/dm-joint-phase-v2/manuscript_dm_catalog.csv in
the Faber2026 parent repo). Reading the yaml here would compare CHIME against
CHIME while labeling one side DSA. Regeneration therefore REQUIRES
``--dsa-catalog CSV`` pointing at that catalog, whose ``dsa_dm``/``dsa_sigma``
columns preserve the independent DSA phase-coherence measurement; without the
flag this script raises. The committed crossmatching/dm_provenance.csv predates
the yaml change and remains the arrival-regression-era artifact.

Trust context (CONTEXT.md, Wave 3): DM_obs is revoked until this per-telescope
provenance exists. This table does NOT by itself re-certify DM_obs -- it exposes
exactly what is and is not documented so the owner can decide. Two structural
findings are baked into the data and surfaced here rather than smoothed over:

  1. The DSA-side dm_err is a uniform 0.1 pc/cm^3 placeholder floor
     (configs/bursts.yaml), not a measured uncertainty. delta_dm_sigma uses the
     same sigma_eff = max(quadrature errors, 1.0 pc/cm^3) convention as
     crossmatching/association.py::dm_agreement -- the 1 pc/cm^3 physical floor
     that keeps a sub-pc offset from reading as many-sigma. So delta_dm_sigma
     matches the committed association_report.json ``n_sigma`` (up to sign), and
     a small value means "consistent within the physical DM scale". See ``note``.
  2. Four bursts (whitney, oran, johndoeII, mahi) are CHIME-unconstrained
     (0-2 sub-bands above S/N 4); their CHIME DM is absent and the agreement is
     undefined -- documented-null (the pipeline ran and reported non-detection),
     not UNDOCUMENTED.

Run from a FLITS checkout at the Faber2026 pin:
    conda run -n flits python scripts/build_dm_provenance.py \
        --dsa-catalog ../analysis/dm-joint-phase-v2/manuscript_dm_catalog.csv
"""
from __future__ import annotations

import csv
import json
import pathlib
from math import hypot

ROOT = pathlib.Path(__file__).parents[1]
CHIME_INPUTS = ROOT / "crossmatching/chime_side_inputs.json"
OUT = ROOT / "crossmatching/dm_provenance.csv"

# Physical DM tolerance floor on the combined sigma, matching
# crossmatching/association.py::dm_agreement (dm_floor). The CHIME arrival
# regression returns a statistical sigma far below the ~1 pc/cm^3 scale at which
# a DM difference is physically meaningful, so without the floor a sub-pc offset
# reads as many-sigma. sigma_eff = max(quadrature errors, DM_FLOOR).
DM_FLOOR = 1.0

# DSA-side DM must come from the manuscript DM catalog's dsa columns (see the
# SOURCE CHANGE note in the module docstring); bursts.yaml is no longer an
# independent DSA source.
DM_DSA_METHOD = "DSA-110 phase-coherence measurement (manuscript DM catalog)"
DM_DSA_SOURCE = (
    "manuscript_dm_catalog.csv [dsa_dm, dsa_sigma] via --dsa-catalog; "
    "bursts.yaml mirrors the adopted CHIME-primary values since 2026-07-27 "
    "and is not an independent DSA source"
)
# CHIME extraction artifacts (chime_dm_final.json, grid NPZ) live off-repo on
# h17 and are not pinned (DM audit gap #2); the pinned repo has only the summary
# JSON + audit note. Cite that chain honestly.
CHIME_SOURCE = (
    "crossmatching/chime_side_inputs.json (audit: .agents/audit-chime-side-dm.md); "
    "extraction artifacts chime_dm_final.json + grid NPZ off-repo on "
    "h17:/data/... (not pinned)"
)

FIELDS = [
    "nickname",
    "dm_dsa", "dm_dsa_err", "dm_dsa_method", "dm_dsa_source",
    "dm_chime", "dm_chime_err", "dm_chime_method", "dm_chime_source",
    "delta_dm", "delta_dm_sigma",
    "dm_confidence", "note",
]


def _load_dsa_catalog(path: pathlib.Path) -> dict[str, tuple[float, float]]:
    # Case-insensitive nick keys (catalog uses johndoeII, other sources johndoeii).
    with open(path, newline="") as fh:
        return {r["nick"].lower(): (float(r["dsa_dm"]), float(r["dsa_sigma"]))
                for r in csv.DictReader(fh)}


def build_rows(dsa_catalog: pathlib.Path | None = None) -> list[dict]:
    if dsa_catalog is None:
        raise RuntimeError(
            "configs/bursts.yaml no longer carries the independent DSA DM "
            "measurement: since 2026-07-27 it mirrors the adopted CHIME-primary "
            "catalog (analysis/dm-joint-phase-v2/manuscript_dm_catalog.csv). "
            "Regenerating dm_provenance.csv from the yaml would compare CHIME "
            "against CHIME while labeling one side DSA. Pass --dsa-catalog "
            "pointing at that catalog to source dsa_dm/dsa_sigma correctly."
        )
    chime = json.load(CHIME_INPUTS.open())
    dsa = _load_dsa_catalog(dsa_catalog)

    rows: list[dict] = []
    for rec in chime:  # preserves the sample (chronological-by-MJD) order
        nick = rec["name"]
        dm_dsa, dm_dsa_err = dsa[nick.lower()]
        dm_chime = rec.get("dm_chime")
        dm_chime_err = rec.get("dm_chime_err")
        method = rec.get("method") or "UNDOCUMENTED"

        if dm_chime is not None and dm_chime_err is not None:
            delta = dm_chime - dm_dsa  # signed: CHIME - DSA (negative = CHIME below)
            # sigma_eff floors the quadrature error at the 1 pc/cm^3 physical
            # scale (association.py convention); |delta_dm_sigma| == report n_sigma.
            sigma_eff = max(hypot(dm_chime_err, dm_dsa_err), DM_FLOOR)
            delta_s = f"{delta:+.4f}"
            sigma_s = f"{delta / sigma_eff:+.3f}"
            note = ("sigma_eff = max(quadrature, 1.0 floor); DSA err is the "
                    "catalog dsa_sigma (phase-coherence measurement)")
        else:
            # CHIME-unconstrained: documented non-detection, agreement undefined.
            delta_s = ""
            sigma_s = ""
            note = "CHIME DM unconstrained (agreement undefined): " + (
                rec.get("dm_status") or "no dm_status recorded"
            )

        rows.append({
            "nickname": nick,
            "dm_dsa": f"{dm_dsa:.4f}",
            "dm_dsa_err": f"{dm_dsa_err:.4f}",
            "dm_dsa_method": DM_DSA_METHOD,
            "dm_dsa_source": DM_DSA_SOURCE,
            "dm_chime": "" if dm_chime is None else f"{dm_chime:.4f}",
            "dm_chime_err": "" if dm_chime_err is None else f"{dm_chime_err:.4f}",
            "dm_chime_method": method,
            "dm_chime_source": CHIME_SOURCE,
            "delta_dm": delta_s,
            "delta_dm_sigma": sigma_s,
            "dm_confidence": rec.get("dm_confidence", ""),
            "note": note,
        })
    return rows


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dsa-catalog", type=pathlib.Path, default=None,
                    help="manuscript_dm_catalog.csv providing independent "
                         "dsa_dm/dsa_sigma (required to regenerate)")
    args = ap.parse_args()
    rows = build_rows(args.dsa_catalog)
    with OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    n_constrained = sum(1 for r in rows if r["delta_dm"])
    print(f"wrote {OUT.relative_to(ROOT)} ({len(rows)} rows; "
          f"{n_constrained} with CHIME-DSA agreement, "
          f"{len(rows) - n_constrained} CHIME-unconstrained)")


if __name__ == "__main__":
    main()
