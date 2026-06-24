#!/usr/bin/env python
"""Manuscript figure: the PBF systematic on the scattering index across the sample.

Dumbbell per burst -- alpha under the (unphysical) per-band mixed PBF -> alpha
under the adopted single-exponential PBF. Well-constrained sightlines barely move
(|d alpha| <= 0.1, dot-sized dumbbells); the poorly-constrained / prior-railed
sightlines are flagged. Values are the all-exp ladder campaign results
(joint_ladder/{LADDER_SUMMARY,ALLEXP_PBF_RUN}.md, 2026-06-24); zach is its
single-component value (its C2D3 was rejected by the fixed-s2 grid).

  python fig_alpha_pbf_systematic.py     # writes alpha_pbf_systematic.{svg,pdf}
"""

import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = os.environ.get("FIG_OUT", ".")

# (TNS, alpha_mixed, alpha_exp, err_minus, err_plus, class)
#   robust   = well-constrained single-screen, |d alpha| <= 0.1
#   marginal = wide posterior reaching a bound
#   railed   = hard-pinned at the alpha=1.5 prior floor under both PBFs
DATA = [
    ("FRB 20230913A", 1.50, 1.504, 0.05, 0.05, "railed"),
    ("FRB 20230814B", 1.58, 1.573, 0.05, 0.05, "railed"),
    ("FRB 20240229A", 2.40, 2.396, 0.04, 0.04, "robust"),
    ("FRB 20221203A", 2.56, 2.558, 0.04, 0.04, "robust"),
    ("FRB 20220506D", 2.69, 2.662, 0.16, 0.16, "robust"),
    ("FRB 20240203A", 3.28, 3.286, 0.04, 0.04, "robust"),
    ("FRB 20220207C", 3.32, 3.319, 0.013, 0.013, "robust"),  # zach, single-comp
    ("FRB 20230307A", 3.33, 3.426, 0.05, 0.05, "robust"),
    ("FRB 20240122A", 3.80, 3.17, 1.18, 1.47, "marginal"),  # mahi, wide->floor
    ("FRB 20230325A", 4.36, 4.356, 0.04, 0.04, "robust"),
    ("FRB 20220310F", 5.21, 5.12, 0.17, 0.17, "robust"),
    ("FRB 20221113A", 5.39, 5.48, 1.98, 0.42, "marginal"),  # isha, wide->upper
]
assert len(DATA) == 12
assert (
    max(abs(e - m) for _, m, e, *_, c in [(d[0], d[1], d[2], d[5]) for d in DATA] if c == "robust")
    <= 0.1
)

COL = {"robust": "#1b7837", "marginal": "#e08214", "railed": "#888888"}
LBL = {
    "robust": r"well-constrained ($|\Delta\alpha|\leq0.1$)",
    "marginal": "marginal (wide posterior)",
    "railed": r"prior-railed ($\alpha\to1.5$)",
}

rows = sorted(DATA, key=lambda r: r[2])  # by alpha_exp
y = np.arange(len(rows))

fig, ax = plt.subplots(figsize=(7.2, 6.4))
ax.axvspan(1.5 - 0.02, 1.7, color="0.92", zorder=0)  # near the lower prior rail
ax.axvline(4.4, ls="--", lw=1.2, color="0.45", zorder=1)
ax.text(4.4, len(rows) - 0.3, " Kolmogorov", fontsize=9, color="0.4", va="top")

for i, (tns, am, ae, em, ep, cls) in enumerate(rows):
    c = COL[cls]
    ax.plot([am, ae], [i, i], "-", color=c, lw=1.5, alpha=0.6, zorder=2)
    ax.plot(am, i, "o", mfc="white", mec=c, mew=1.4, ms=7, zorder=3)  # mixed PBF
    ax.errorbar(
        ae, i, xerr=[[em], [ep]], fmt="o", color=c, ms=7, capsize=2.5, elinewidth=1.1, zorder=4
    )  # adopted all-exp PBF

ax.set_yticks(y)
ax.set_yticklabels([r[0] for r in rows], fontsize=9)
ax.set_xlabel(r"scattering index $\alpha$")
ax.set_xlim(1.2, 6.1)
ax.set_ylim(-0.6, len(rows) - 0.4)
ax.set_title(
    "PBF systematic on the scattering index: per-band (open) "
    r"$\to$ single-exponential (filled)",
    fontsize=10.5,
)

handles = [
    Line2D(
        [], [], marker="o", mfc="white", mec="0.3", mew=1.4, ls="", label="per-band PBF (mixed)"
    ),
    Line2D([], [], marker="o", color="0.3", ls="", label="single-exponential PBF (adopted)"),
]
handles += [
    Line2D([], [], marker="o", color=COL[k], ls="", label=LBL[k])
    for k in ("robust", "marginal", "railed")
]
ax.legend(handles=handles, fontsize=8, loc="lower right", framealpha=0.9)
fig.tight_layout()
for ft in ("svg", "pdf"):
    fp = f"{OUT}/alpha_pbf_systematic.{ft}"
    fig.savefig(fp, bbox_inches="tight")
    print("wrote", fp)
