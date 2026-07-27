#!/usr/bin/env python
"""Final joint-fit campaign table (amended 2026-07-17): per burst, per instrument.

Columns (owner spec): chosen component counts, TOA +/- err, alpha, beta, tau_1ghz,
delta_dm per band, resolution factors (f/t) + achieved peak S/N, residual max per
band, and OLD-vs-NEW where an OLD fit exists.

TOA REFERENCE (uniform, stated): the FLUENCE-WEIGHTED CENTROID of the component
arrival times per band, weights = OLS-recovered per-component spectral fluence at
the posterior median (from the jointmodel dump), error = posterior spread of the
centroid (component t0's varied over the equal-weight posterior, weights fixed).
Reduces exactly to t0 +/- err for a single-component band. Component-count changes
vs OLD, and bursts where CHIME and DSA resolve DIFFERENT counts (matched reference
most delicate), are flagged explicitly.

  FLITS_RUNS=~/Developer/scratch/flits-local-runs \
  CAMPAIGN_LOGS=<...>/campaign_A1 python build_toa_table.py
"""
from __future__ import annotations

import csv
import json
import os
import re
from pathlib import Path

import numpy as np
from dynesty.utils import resample_equal
from scat_analysis.turbulence import alpha_from_beta

RUNS = Path(os.environ.get("FLITS_RUNS", os.path.expanduser("~/Developer/scratch/flits-local-runs")))
JOINT = RUNS / "data/joint"
OLD = Path(os.environ.get("OLD_FITS", "/private/tmp/claude-501/-Users-jakobfaber-Developer-repos-github-com-jakobtfaber-Faber2026/a573656f-9ea8-4ebd-aab6-58332c63c659/scratchpad/fits_OLD_campaign"))
LOGS = Path(os.environ.get("CAMPAIGN_LOGS", "/private/tmp/claude-501/-Users-jakobfaber-Developer-repos-github-com-jakobtfaber-Faber2026/a573656f-9ea8-4ebd-aab6-58332c63c659/scratchpad/campaign_A1"))
HERE = Path(__file__).resolve().parent

# burst -> (new_tag, old_tag). Chosen counts = joint_ladder/_figs.py `chosen` map,
# beta-native re-fit. whitney corrected base->C2D2 per owner. Component-count deltas
# vs OLD (hamilton sharedzeta->C4D1, zach C1D1->C2D3) are surfaced as flags below.
BURSTS = [
    ("freya", "", "_sharedzeta"),
    ("casey", "", "_sharedzeta"),
    ("chromatica", "", "_sharedzeta"),
    ("wilhelm", "", "_sharedzeta"),
    ("hamilton", "_C4D1", "_sharedzeta"),
    ("mahi", "_C1D1", "_C1D1"),
    ("oran", "_C2D1", "_C2D1"),
    ("isha", "_C2D1", "_C2D1"),
    ("whitney_fine", "_C2D2", "_C2D2"),
    ("johndoeII", "_C2D2", "_C2D2"),
    ("phineas", "_C3D3", "_C3D3"),
    ("zach", "_C2D3", "_C1D1"),
]

# Bands excluded from the TOA analysis for data-availability reasons (not physics):
# isha's DSA product is unusable for timing (owner decision), so its DSA-band TOA
# is not reported.
EXCLUDE_D = {"isha"}

# beta_native alpha_from_beta clamps to exactly 4.0 for beta >= 4 - BETA_EXP_EPS
# (= 3.98): inside that exponential-PBF regime alpha is pinned at 4.0, so the
# derived alpha is NOT a free measurement. The one-sided limits below are the
# honest statement of the constraint; the UNCLAMPED alpha constraint comes from
# the relaxed-alpha A/B, not this table.
#
# Two rails exist. HIGH: beta -> 4 (Gaussian/thin-screen ceiling), alpha -> 4 --
# report beta > lo95, alpha < hi95. LOW: beta -> 3 (Kolmogorov steep-index floor),
# alpha -> 6 -- the COMPLEMENTARY limit, report beta < hi95, alpha > lo95 (phineas).
BETA_RAIL_HI = 3.98   # >= this fraction of samples above -> rails at beta=4
BETA_RAIL_LO = 3.05   # <= this fraction of samples below -> rails at beta=3


def rail_limits(samp_npz: Path):
    """Corner-aware one-sided 95% limits on beta and (per-sample derived) alpha.

    Returns None if the samples are absent. ``corner`` is "high" (beta=4 ceiling),
    "low" (beta=3 floor / alpha=6), or "resolved". ``limit_str`` is the manuscript
    phrasing for the relevant tail; a two-sided posterior reports the median with
    its 5-95 spread instead of a spurious rail.
    """
    if not samp_npz.exists():
        return None
    z = np.load(samp_npz, allow_pickle=True)
    names = list(z["param_names"])
    if "beta" not in names:
        return None
    bi = names.index("beta")
    eq = resample_equal(z["samples"], z["weights"])
    b = eq[:, bi]
    a = np.array([alpha_from_beta(float(min(max(bb, 2.001), 4.0))) for bb in b])
    b_lo95, b_hi95 = float(np.percentile(b, 5)), float(np.percentile(b, 95))
    a_lo95, a_hi95 = float(np.percentile(a, 5)), float(np.percentile(a, 95))
    frac_hi = float(np.mean(b > BETA_RAIL_HI))
    frac_lo = float(np.mean(b < BETA_RAIL_LO))
    if frac_hi >= 0.90:
        corner, limit_str = "high", f"beta > {b_lo95:.3f} (95%), alpha < {a_hi95:.3f} (95%)"
    elif frac_lo >= 0.90:
        corner, limit_str = "low", f"beta < {b_hi95:.3f} (95%), alpha > {a_lo95:.3f} (95%)"
    else:
        corner, limit_str = "resolved", (
            f"beta = {float(np.median(b)):.3f} [{b_lo95:.3f}, {b_hi95:.3f}], "
            f"alpha = {float(np.median(a)):.3f} [{a_lo95:.3f}, {a_hi95:.3f}]"
        )
    return dict(
        beta_med=float(np.median(b)), alpha_med=float(np.median(a)),
        beta_lo95=b_lo95, beta_hi95=b_hi95, alpha_lo95=a_lo95, alpha_hi95=a_hi95,
        frac_beta_hi=frac_hi, frac_beta_lo=frac_lo,
        corner=corner, railed=(corner != "resolved"), limit_str=limit_str,
    )


def load(fp: Path):
    return json.load(open(fp)) if fp.exists() else None


def pget(d, k):
    if not d:
        return None
    p = d.get("percentiles", {}).get(k)
    if p:
        return p["median"], p.get("err_minus", 0.0), p.get("err_plus", 0.0)
    if k in d and isinstance(d[k], dict) and "median" in d[k]:
        return d[k]["median"], d[k].get("err_minus", 0.0), d[k].get("err_plus", 0.0)
    return None


def t0_cols(names, band):
    """Column indices for a band's t0 params (t0_C or t0_C1..t0_C{n}), in order."""
    idx = [(n, i) for i, n in enumerate(names) if re.fullmatch(rf"t0_{band}\d*", n)]
    idx.sort(key=lambda x: (len(x[0]), x[0]))  # t0_C before t0_C1..; numeric order otherwise
    return [i for _, i in idx], [n for n, _ in idx]


def centroid_toa(samples_npz: Path, fluence, band):
    """Fluence-weighted centroid TOA (median, err_minus, err_plus) + per-component
    t0 medians, from the equal-weight posterior. Weights fixed at median fluence."""
    if not samples_npz.exists():
        return None, []
    z = np.load(samples_npz, allow_pickle=True)
    names = list(z["param_names"])
    cols, cnames = t0_cols(names, band)
    if not cols:
        return None, []
    eq = resample_equal(z["samples"], z["weights"])   # (E, P) equal-weight
    t0s = eq[:, cols]                                   # (E, k)
    w = np.ones(len(cols), float)
    if fluence is not None and len(fluence) == len(cols):
        w = np.clip(np.asarray(fluence, float), 0.0, None)
        if w.sum() <= 0:
            w = np.ones(len(cols), float)
    w = w / w.sum()
    cen = t0s @ w                                       # (E,)
    med, lo, hi = np.percentile(cen, [50, 16, 84])
    comp_meds = [(cnames[j], float(np.median(t0s[:, j])), float(w[j])) for j in range(len(cols))]
    return (float(med), float(med - lo), float(hi - med)), comp_meds


def caption_from_log(burst):
    fp = LOGS / f"{burst}.log"
    if not fp.exists():
        return {}
    txt = fp.read_text()
    out = {}
    for band in ("CHIME", "DSA"):
        m = re.search(rf"AUTO-TF {band}\s*:\s*(.*)", txt)
        if m:
            c = m.group(1).strip()
            g = dict(caption=c)
            for key, pat in (("ff", r"f(\d+)/"), ("tf", r"/t(\d+)\)"),
                             ("win_ms", r"window ([\d.]+) ms"), ("snr", r"peak S/N (\d+)/px"),
                             ("nch", r"(\d+) ch"), ("dt_us", r"dt=([\d.]+) us")):
                mm = re.search(pat, c)
                if mm:
                    g[key] = mm.group(1)
            out[band] = g
    return out


def resid_from_json(burst, tag):
    fp = JOINT / f"{burst}_jointmodel{tag}_resid.json"
    if not fp.exists():
        return None
    return json.load(open(fp))


def dump_arrays(burst, tag):
    npz = JOINT / f"{burst}_jointmodel{tag}.npz"
    if not npz.exists():
        return None
    z = np.load(npz, allow_pickle=True)
    return dict(
        chi2C=float(z["chi2C"]) if "chi2C" in z else None,
        chi2D=float(z["chi2D"]) if "chi2D" in z else None,
        fluenceC=z["fluenceC"] if "fluenceC" in z else None,
        fluenceD=z["fluenceD"] if "fluenceD" in z else None,
        nC=int(z["nC"]) if "nC" in z else 1,
        nD=int(z["nD"]) if "nD" in z else 1,
    )


def fmt(v, u=1e3):
    return "--" if v is None else f"{v[0]:+.4f} (+{v[2]:.4f}/-{v[1]:.4f})"


def shape_inflate(toa, chi2, resid):
    """Inflate a band's TOA (+/-)error by sqrt(chi2_red) when its residual is a
    shape mismatch. A single-component model that is shape-misspecified has
    chi2_red > 1, and the posterior-spread error understates the true timing
    uncertainty by ~sqrt(chi2_red) (the standard EFAC correction). Returns the
    TOA unchanged for well-fit bands."""
    if toa is None or not resid or not resid.get("shape_mismatch") or not chi2:
        return toa, False
    s = float(np.sqrt(max(float(chi2), 1.0)))
    return (toa[0], toa[1] * s, toa[2] * s), (s > 1.0)


def main():
    rows = []
    print("=" * 120)
    print("JOINT-FIT CAMPAIGN (amended): NEW = beta-native [3,4], AUTO S/N-driven prep, evidence-selected counts")
    print("TOA = fluence-weighted centroid of component arrival times per band (uniform reference)")
    print("=" * 120)
    for burst, ntag, otag in BURSTS:
        new = load(JOINT / f"{burst}_joint_fit{ntag}.json")
        old = load(OLD / f"{burst}_joint_fit{otag}.json")
        cap = caption_from_log(burst)
        da = dump_arrays(burst, ntag)
        rj = resid_from_json(burst, ntag)
        nC = (new or {}).get("components_C", 1)
        nD = (new or {}).get("components_D", 1)
        missing = "" if new else "  [NEW FIT MISSING]"
        # TOA centroids (NEW)
        samp = JOINT / f"{burst}_joint_samples{ntag}.npz"
        toaC, compC = centroid_toa(samp, (da or {}).get("fluenceC"), "C")
        toaD, compD = centroid_toa(samp, (da or {}).get("fluenceD"), "D")
        rl = rail_limits(samp)

        flags = []
        # DSA band excluded from timing for data-availability reasons.
        if burst in EXCLUDE_D:
            toaD, compD = None, []
            flags.append("DSA TOA EXCLUDED (data availability, not physics)")
        # Inflate shape-mismatch band TOA errors by sqrt(chi2_red) (EFAC).
        toaC, infC = shape_inflate(toaC, (da or {}).get("chi2C"), (rj or {}).get("C"))
        toaD, infD = shape_inflate(toaD, (da or {}).get("chi2D"), (rj or {}).get("D"))
        if infC:
            flags.append("CHIME TOA error INFLATED x sqrt(chi2_red) (shape mismatch)")
        if infD:
            flags.append("DSA TOA error INFLATED x sqrt(chi2_red) (shape mismatch)")
        if rl and rl["railed"]:
            ceiling = "beta=4 ceiling (alpha=4)" if rl["corner"] == "high" else "beta=3 floor (alpha=6)"
            quote = "alpha=4.0" if rl["corner"] == "high" else "alpha=6.0"
            flags.append(
                f"beta AT {ceiling}: report LIMIT {rl['limit_str']} — do NOT quote {quote} as a measurement"
            )
        if nC != nD:
            flags.append(f"CHIME/DSA resolve DIFFERENT counts (C{nC} vs D{nD}) — matched-ref delicate")
        oc = (old or {}).get("components_C", 1)
        od = (old or {}).get("components_D", 1)
        if old and (oc != nC or od != nD):
            flags.append(f"component count CHANGED vs OLD (C{oc}D{od}->C{nC}D{nD}) — TOA may shift materially")
        for band, r in (("CHIME", (rj or {}).get("C")), ("DSA", (rj or {}).get("D"))):
            if r and r.get("escalate"):
                flags.append(f"{band} residual ESCALATE (resid_max {r['resid_prof_max']:+.1f}s, "
                             f"{r['n_contig_5sig']}bin) — ignored component; refit at higher count")
            elif r and r.get("shape_mismatch"):
                flags.append(f"{band} residual SHAPE-MISMATCH (+/-{r['resid_prof_max_abs']:.1f}s dipole) "
                             f"— bright-pulse shape/resolution, NOT a missing component")

        print(f"\n----- {burst}  NEW C{nC}D{nD}{missing} -----")
        for band, g in (("CHIME", cap.get("CHIME")), ("DSA", cap.get("DSA"))):
            if g:
                print(f"   {band:5s} resolution: f{g.get('ff','?')}/t{g.get('tf','?')}  "
                      f"{g.get('nch','?')}ch  dt={g.get('dt_us','?')}us  win={g.get('win_ms','?')}ms  "
                      f"peakS/N={g.get('snr','?')}/px")
        print(f"   chi2_red  CHIME={(da or {}).get('chi2C')}  DSA={(da or {}).get('chi2D')}")
        if rj:
            rc, rd = rj.get("C", {}), rj.get("D", {})
            print(f"   residual_max  CHIME={rc.get('resid_prof_max'):+.2f}s (contig5s={rc.get('n_contig_5sig')})  "
                  f"DSA={rd.get('resid_prof_max'):+.2f}s (contig5s={rd.get('n_contig_5sig')})")
        for key, lab in [("tau_1ghz", "tau_1GHz(ms)"), ("alpha", "alpha"), ("beta", "beta"),
                         ("delta_dm_C", "dDM_C(pc/cc)"), ("delta_dm_D", "dDM_D(pc/cc)")]:
            print(f"   {lab:14s} OLD {fmt(pget(old,key)):32s}  NEW {fmt(pget(new,key))}")
        if rl:
            tag_r = {"high": "RAILED beta=4", "low": "RAILED beta=3", "resolved": "resolved"}[rl["corner"]]
            print(f"   scatter {rl['limit_str']}"
                  f"   [{tag_r}; frac_beta>3.98={rl['frac_beta_hi']:.2f}, frac_beta<3.05={rl['frac_beta_lo']:.2f}]")
        print(f"   TOA_CHIME (centroid, ms)  NEW {fmt(toaC)}")
        for nm, m, w in compC:
            print(f"       component {nm:7s} t0={m:+.4f} ms  (fluence weight {w:.2f})")
        print(f"   TOA_DSA   (centroid, ms)  NEW {fmt(toaD)}")
        for nm, m, w in compD:
            print(f"       component {nm:7s} t0={m:+.4f} ms  (fluence weight {w:.2f})")
        if flags:
            for fl in flags:
                print(f"   >> FLAG: {fl}")

        def snr(band):
            g = cap.get(band) or {}
            return g.get("snr")
        rcC = (rj or {}).get("C", {}) or {}
        rcD = (rj or {}).get("D", {}) or {}
        tau_n = pget(new, "tau_1ghz")
        al_n = pget(new, "alpha")
        rows.append(dict(
            burst=burst, comp=f"C{nC}D{nD}",
            chime_ff_tf=f"f{(cap.get('CHIME') or {}).get('ff','?')}/t{(cap.get('CHIME') or {}).get('tf','?')}",
            dsa_ff_tf=f"f{(cap.get('DSA') or {}).get('ff','?')}/t{(cap.get('DSA') or {}).get('tf','?')}",
            chime_peaksnr=snr("CHIME"), dsa_peaksnr=snr("DSA"),
            chi2_C=(da or {}).get("chi2C"), chi2_D=(da or {}).get("chi2D"),
            tau_ms=(tau_n[0] if tau_n else None),
            alpha=(al_n[0] if al_n else None),
            beta_med=(rl["beta_med"] if rl else None),
            beta_lo95=(rl["beta_lo95"] if rl else None),
            beta_hi95=(rl["beta_hi95"] if rl else None),
            alpha_lo95=(rl["alpha_lo95"] if rl else None),
            alpha_hi95=(rl["alpha_hi95"] if rl else None),
            rail_corner=(rl["corner"] if rl else None),
            railed=(rl["railed"] if rl else None),
            scatter_limit=(rl["limit_str"] if rl else None),
            dsa_excluded=(burst in EXCLUDE_D),
            toa_C_ms=(toaC[0] if toaC else None),
            toa_C_err_ms=((toaC[1] + toaC[2]) / 2 if toaC else None),
            toa_D_ms=(toaD[0] if toaD else None),
            toa_D_err_ms=((toaD[1] + toaD[2]) / 2 if toaD else None),
            resid_max_C=rcC.get("resid_prof_max"), resid_max_D=rcD.get("resid_prof_max"),
            escalate_C=rcC.get("escalate"), escalate_D=rcD.get("escalate"),
            flags="; ".join(flags),
        ))
    csv_fp = HERE / "joint_tf_toa_table.csv"
    with open(csv_fp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {csv_fp}")


if __name__ == "__main__":
    main()
