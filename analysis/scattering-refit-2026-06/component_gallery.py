#!/usr/bin/env python
"""Per-burst component-screening gallery: CHIME | DSA, each a dedispersed
waterfall (freq x time, on-pulse cropped) over its frequency-collapsed profile.
No model overlay -- this is for COUNTING temporal sub-components by eye before
deciding C1 vs C2 fits. Reuses the exact data prep as run_joint_fit/joint_ppc.

  python component_gallery.py <burst>
"""

import os
import sys

import matplotlib
import numpy as np
import yaml

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.environ.get("FLITS_REPO", "/home/jfaber/flits/dsa110-FLITS")
RUNS = os.environ.get("FLITS_RUNS", "/central/scratch/jfaber/flits-runs")
sys.path.insert(0, f"{REPO}/scattering")
from scat_analysis.config_utils import load_telescope_block
from scat_analysis.pipeline.io import BurstDataset


def prepare(cfg_path, name):
    cfg = yaml.safe_load(open(cfg_path))
    tel = load_telescope_block(cfg["telcfg_path"], cfg["telescope"])
    ds = BurstDataset(
        cfg["path"],
        f"{RUNS}/data/joint",
        name=name,
        telescope=tel,
        f_factor=int(cfg["f_factor"]),
        t_factor=int(cfg["t_factor"]),
        outer_trim=float(cfg.get("outer_trim", 0.15)),
        onpulse_crop=os.environ.get("FLITS_ONPULSE_CROP", "1") == "1",
        onpulse_pad_factor=float(os.environ.get("FLITS_ONPULSE_PAD", "0.5")),
    )
    return ds.model


def panel(axw, axp, m, title):
    d = np.asarray(m.data, float)  # (freq, time), freq ascending
    t = np.asarray(m.time, float)
    # downsample freq for display only if very tall
    step = max(1, d.shape[0] // 512)
    disp = d[::step]
    vmax = np.nanpercentile(disp, 99.5)
    axw.imshow(
        disp,
        aspect="auto",
        origin="lower",
        extent=[t[0], t[-1], 0, d.shape[0]],
        cmap="viridis",
        vmin=np.nanpercentile(disp, 5),
        vmax=vmax,
    )
    axw.set_title(title)
    axw.set_ylabel("freq chan")
    prof = np.nansum(d, axis=0)
    axp.plot(t, prof, "k", lw=0.9)
    axp.axhline(0, color="0.6", lw=0.5)
    axp.set_xlabel("time (ms)")
    axp.set_ylabel("flux (a.u.)")
    axp.margins(x=0)


def main():
    b = sys.argv[1]
    mC = prepare(f"{RUNS}/configs/{b}_chime_run.yaml", f"{b}_chime")
    mD = prepare(f"{RUNS}/configs/{b}_dsa_run.yaml", f"{b}_dsa")
    fig, ax = plt.subplots(2, 2, figsize=(13, 6), gridspec_kw={"height_ratios": [2, 1]})
    panel(ax[0, 0], ax[1, 0], mC, f"{b}  CHIME (400-800 MHz)")
    panel(ax[0, 1], ax[1, 1], mD, f"{b}  DSA (1.31-1.50 GHz)")
    fig.suptitle(f"{b}: component screen -- count temporal sub-pulses per band", y=1.0)
    fig.tight_layout()
    fp = f"{RUNS}/data/joint/{b}_components.png"
    fig.savefig(fp, dpi=120, bbox_inches="tight")
    print(f"wrote {fp}")


if __name__ == "__main__":
    main()
