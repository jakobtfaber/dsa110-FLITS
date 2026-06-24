#!/usr/bin/env python
"""Manuscript figure: whitney (FRB 20220310F) DSA-band multiplicity contrast.

Two rows -- single DSA component (rails alpha->1.5 because one scattered pulse
cannot cover both 0.34 ms sub-pulses) vs two DSA components (unrails to
alpha=5.1, white residual). DSA band only: that is where the multiplicity lives.
Reconstruction (gain-marginal per-channel solve) is the same as model_2d.py.

  python fig_whitney_multiplicity.py
Outputs whitney_multiplicity.{svg,pdf} into the overleaf figures dir.
"""

import json
import os
import sys
from dataclasses import replace

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.environ.get("FLITS_REPO")
RUNS = os.environ.get("FLITS_RUNS")
OUT = os.environ.get("FIG_OUT", ".")
sys.path.insert(0, RUNS)
sys.path.insert(0, f"{REPO}/scattering")
from run_joint_fit import prepare
from scat_analysis.burstfit import FRBParams
from scat_analysis.burstfit_joint import _gain_marginal_multi_band

BURST = "whitney_fine"
# (tag in JSON filename, n DSA components, row label) -- both all-exponential PBF.
ROWS = [
    ("C2D1_s2-1_pbf-exp-exp", 1, r"1 DSA component: $\alpha\!\to\!1.5$ (railed)"),
    ("C2D2_pbf-exp-exp", 2, r"2 DSA components: $\alpha=5.1$"),
]


def load(tag):
    J = json.load(open(f"{RUNS}/data/joint/{BURST}_joint_fit_{tag}.json"))
    P = {k: v["median"] for k, v in J["percentiles"].items()}
    a = J["alpha"]["median"]
    return P, J["percentiles"]["tau_1ghz"]["median"], a


def band_params(P, n, tau, alpha, ddm):
    return [
        FRBParams(
            c0=1.0,
            t0=P[f"t0_D{i}"],
            gamma=0.0,
            zeta=P[f"zeta_D{i}"],
            tau_1ghz=tau,
            alpha=alpha,
            delta_dm=ddm,
        )
        for i in range(1, n + 1)
    ]


def solve(model, ps):
    valid = model.valid
    Ks = np.stack([model(replace(p, c0=1.0, gamma=0.0), "M3", freq_subset=valid) for p in ps])
    d = model.data[valid]
    sig = np.clip(model.noise_std[valid], 1e-9, None)
    var = sig**2
    b = np.einsum("nft,ft->fn", Ks, d)
    M = np.einsum("nft,mft->fnm", Ks, Ks)
    _, diag = _gain_marginal_multi_band(model, ps, ["M3"] * len(ps), s2=None)
    A = M + (var / diag["s2"])[:, None, None] * np.eye(Ks.shape[0])[None]
    g = np.linalg.solve(A, b[:, :, None])[:, :, 0]
    pred = np.einsum("fn,nft->ft", g, Ks)
    resid = (d - pred) / sig[:, None]
    chi2 = float(np.sum(resid**2) / resid.size)
    return d, pred, resid, model.time, model.freq[valid], chi2


model = prepare(f"{RUNS}/configs/{BURST}_dsa_run.yaml", f"{BURST}_dsa", f"{RUNS}/data/joint")[0]

fig, ax = plt.subplots(2, 3, figsize=(11, 5.4), sharex=True, sharey=True)
for r, (tag, nD, label) in enumerate(ROWS):
    P, tau, alpha = load(tag)
    d, pred, resid, t, f, chi2 = solve(
        model, band_params(P, nD, tau, alpha, P.get("delta_dm_D", 0.0))
    )
    ext = [t[0], t[-1], f[0], f[-1]]
    vmax = np.nanpercentile(d, 99.5)
    vmin = np.nanpercentile(d, 2)
    panels = [
        (d, "data", "viridis", dict(vmin=vmin, vmax=vmax)),
        (pred, f"model ({nD} comp)", "viridis", dict(vmin=vmin, vmax=vmax)),
        (resid, "residual", "RdBu_r", dict(vmin=-4, vmax=4)),
    ]
    for c, (img, title, cmap, kw) in enumerate(panels):
        im = ax[r, c].imshow(img, aspect="auto", origin="lower", extent=ext, cmap=cmap, **kw)
        if r == 0:
            ax[r, c].set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax[r, c], pad=0.01, fraction=0.045)
    ax[r, 2].text(
        0.04,
        0.92,
        rf"$\chi^2_\nu={chi2:.2f}$",
        transform=ax[r, 2].transAxes,
        fontsize=10,
        va="top",
        ha="left",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.6", alpha=0.85),
    )
    ax[r, 0].set_ylabel(f"{label}\nfreq (GHz)", fontsize=9)
for c in range(3):
    ax[1, c].set_xlabel("time (ms)")
fig.suptitle(
    "FRB 20220310F (DSA 1.31--1.50 GHz): scattering index rails when the "
    "second sub-pulse is unmodeled",
    fontsize=11,
)
fig.tight_layout()
for ftype in ("svg", "pdf"):
    fp = f"{OUT}/whitney_multiplicity.{ftype}"
    fig.savefig(fp, bbox_inches="tight")
    print("wrote", fp)
