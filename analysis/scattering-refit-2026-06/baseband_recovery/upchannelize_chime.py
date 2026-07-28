#!/usr/bin/env python3
"""CHIME baseband upchannelization for the resolution-limited co-detection sightlines.

Runs INSIDE the `chimefrb/baseband-analysis:latest` docker image on h17 (lxd110h17), which carries
baseband_analysis 1.9.0 + the CADC `vos` client. h17 reaches the CHIME baseband store on CANFAR/arc
directly (verified: `vls arc:projects/chime_frb/...` works in-container with ~/.ssl/cadcproxy.pem),
so there is NO CANFAR Science-Platform / Harbor dependency. Campaign mode uses only explicitly
staged h17 inputs and never fetches data implicitly. Per target this:
  1. verifies the staged singlebeam_<id>.h5 and records its SHA-256,
  2. coherently dedisperses the complex per-channel baseband at the burst DM,
  3. upchannelizes each 0.390625 MHz CHIME coarse channel by the verified per-target factor,
  4. forms Stokes I and writes an immutable factor/software/container-tagged product directory.

WHY coherent dedispersion + PFB upchannelization (not a cheap incoherent rechannel):
  The scintillation measurement is a spectral autocorrelation (ACF) of the time-integrated burst
  spectrum; its diffractive bandwidth Dnu_d is the HWHM of the ACF's central Lorentzian. Two
  systematics counterfeit a scintle and bias Dnu_d if not removed at the baseband level:
    1. Intra-channel dispersive smearing. At these DMs (462-960 pc/cc) the sweep across one CHIME
       coarse channel is many microseconds; upchannelizing a smeared channel imprints a dispersive
       chirp that survives as spurious narrow-band ACF structure. Coherent dedispersion on the raw
       voltages removes the sweep EXACTLY (the only way to recover sub-channel spectral resolution
       without re-smearing). Incoherent (post-detection) dedispersion cannot -- the phase needed to
       de-chirp is gone after squaring.
    2. PFB channel-edge response. The CHIME polyphase-filterbank inverse (baseband_analysis) gives
       the synthesized fine channels a flat passband; a naive FFT-rechannel leaves the coarse-channel
       scallop, a deterministic ripple that contaminates the small-lag ACF exactly where Dnu_d lives.
  The limiting scintle at CHIME is NARROWER than one 0.390625 MHz coarse channel (NE2025/NE2001
  predict sub-channel Dnu_d for all of these), so it is unresolved at native CHIME resolution;
  upchannelizing exposes it -- but only if the channel is de-chirped first.

UPCHAN FACTOR (skeptic-corrected, sized to the DOMINANT/narrower scintle at >=4 ch across its HWHM):
  - casey   U=16  host-dominated  (host Dnu_d 0.187 < MW floor 0.207 MHz) -- the one clean host case.
  - whitney U=16  MW-floor-dominated (resolves the 0.140 MHz floor; native x16 over-resolves, fine).
  - phineas U=16  MW-floor-dominated (0.206 MHz floor; min factor is 8 but native x16 is cleaner).
  - mahi    U=512 MW-floor-dominated (0.0036 MHz floor); only non-smeared because FWHM=24 ms. Uses
            the _upchannel(fftsize=2U) generalization (see _waterfall) -- the SLOW python-loop path.
  - isha    OFF by default: NOT cleanly resolvable -- DSA input gamma railed at the 0.06 MHz fit floor
            AND the dominant scale needs U>=256-512 while the 1.8 ms burst smears to <3 time elements.
            Run only as a lower-confidence upper bound with --run-unresolvable.

API NOTE: baseband_analysis.analysis.waterfall_from_beamformed is BROKEN in this image (v1.9.0): it
feeds upchannel()'s 3-tuple straight into incoherent_dedisp, which does matrix_in.copy() ->
AttributeError on the tuple. So this worker drives the package's own primitives directly for ALL
factors: coherent_dedisp + the internal _upchannel(fftsize=2U, downfreq=2) (which returns
(spec, freq, chan_id); upchan factor U = fftsize/downfreq), forming Stokes I from the complex spectrum
without the incoherent step (coherent dedispersion already de-chirps fully). See _waterfall.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shlex
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

CHIME_COARSE_DF_MHZ = 0.390625  # CHIME coarse channel width (400 MHz / 1024)
CHIME_NATIVE_DT_S = 2.56e-6  # CHIME single-pol baseband sample time
CHIME_NOMINAL_COARSE_CHANNELS = 1024
CHIME_TOP_EDGE_MHZ = 800.1953125
CHIME_BOTTOM_EDGE_MHZ = 400.1953125
SUPPORTED_UPCHANNEL_FACTORS = (16, 32, 64, 128, 256, 512)
INSTRUMENT_ID = "CHIME/FRB"
DM_PROVENANCE_SCHEMA = "faber2026-chime-dm-provenance-v1"
PRODUCT_MANIFEST_SCHEMA = "faber2026-chime-upchannel-product-v1"
PLAN_MANIFEST_SCHEMA = "faber2026-chime-upchannel-campaign-plan-v1"

# All 12 singlebeam .h5 are staged by project ID on h17. Campaign operation
# fails closed when a staged input is absent; it never fetches data implicitly.
LOCAL_H5_DIR = "/data/Faber2026/data/chime-frb"
DEFAULT_OUT_DIR = "/data/research/astrophysics/frbs/chime-dsa-codetections/upchan_codetections"

# id/dm/fwhm_ms from crossmatching/notebook_reproduction_fixture.json (DM also in configs/bursts.yaml).
# Batch 2 (2026-07-07, remaining 6 sightlines): U sized by the same rule -- next power of 2 with
# >=4 fine channels across the NE2025 MW-floor HWHM at 600.19 MHz (U >= 1562.5/floor_kHz), floors
# re-derived with scintillation/ne2025/query_ne2025_scint.py (reproduces the batch-1 floors exactly).
# METHOD EXTENSION vs batch 1: hamilton and johndoeII fail the >=3-time-element gate at their
# frequency-required U, but their time-bandwidth volume W*Dnu_d ~ 10-12 >> 1, so the scintle
# structure IS present in a single upchannelized FFT block ("single-block mode"): the burst lands
# in ~1 output time sample; downstream must take the one on-burst spectrum directly (no
# matched-filter profile gain) and verify the burst does not straddle a block boundary from
# time0 metadata. They stay recoverable:False so the gate is an explicit --run-unresolvable choice.
TARGETS = {
    "casey": {
        "id": "362593221",
        "dm": 491.207,
        "fwhm_ms": 0.1798,
        "upchan": 16,
        "recoverable": True,
        "h5_relpath": "2024/02/29/astro_362593221/singlebeam_362593221.h5",
        "note": "host-dominated (host 0.187 < MW floor 0.207 MHz); cleanest DSA input -- the clean host recovery.",
    },
    "whitney": {
        "id": "215063905",
        "dm": 462.174,
        "fwhm_ms": 0.4865,
        "upchan": 16,
        "recoverable": True,
        "h5_relpath": "2022/03/10/astro_215063905/singlebeam_215063905.h5",
        "note": "MW-floor-dominated (0.140 MHz floor); native x16 (24 kHz) resolves it. Galactic, not host.",
    },
    "phineas": {
        "id": "274819243",
        "dm": 610.274,
        "fwhm_ms": 2.9886,
        "upchan": 16,
        "recoverable": True,
        "h5_relpath": "2023/03/07/astro_274819243/singlebeam_274819243.h5",
        "note": "MW-floor-dominated (0.206 MHz floor); native x16. Long FWHM, time-res ample. Galactic, not host.",
    },
    "mahi": {
        "id": "354049284",
        "dm": 960.128,
        "fwhm_ms": 24.286,
        "upchan": 512,
        "recoverable": True,
        "h5_relpath": "2024/01/22/astro_354049284/singlebeam_354049284.h5",
        "note": "MW-floor-dominated (0.0036 MHz floor) -> U=512 via _upchannel(fftsize=1024); slow. Safe only b/c FWHM=24 ms.",
    },
    "freya": {
        "id": "278720455",
        "dm": 912.4,
        "fwhm_ms": 0.40065937342767416,
        "upchan": 64,
        "recoverable": True,
        "h5_relpath": "2023/03/25/astro_278720455/singlebeam_278720455.h5",
        "note": "MW-floor-sized: NE2025 floor 1.6421 MHz @1.405 GHz -> 38.8 kHz @600 MHz (nu^-4.4, census convention); U=64 (df 6.10 kHz) ~6.4 ch across it. Time: dt=0.328 ms vs intrinsic FWHM 0.401 ms (1.2 el) BUT the CHIME profile is scattering-broadened (tau_600~1.07 ms per the beta co-model production fit) -> >=3 elements across the observed profile. The contested ~0.51 MHz DSA sub-floor candidate (~10 kHz @600) is NOT resolvable at any U inside the time wall (needs U~256, dt 1.31 ms).",
    },
    "zach": {
        "id": "210456524",
        "dm": 262.368,
        "fwhm_ms": 0.964,
        "upchan": 64,
        "recoverable": True,
        "h5_relpath": "2022/02/07/astro_210456524/run_pre_apr2025/singlebeam_210456524.h5",
        "note": "MW-floor-dominated (36.0 kHz floor); U=64 (6.10 kHz) ~5.9 ch. Time: dt=0.328 ms vs "
        "scattering-broadened profile ~2.1 ms (tau600~2.1 from joint fit, alpha=3.66) -> ~6 el.",
    },
    "chromatica": {
        "id": "356959136",
        "dm": 272.664,
        "fwhm_ms": 0.821,
        "upchan": 64,
        "recoverable": True,
        "h5_relpath": "2024/02/03/astro_356959136/old_processed_files/singlebeam_356959136.h5",
        "note": "MW-floor-dominated (35.3 kHz floor); U=64 ~5.8 ch. Time: profile ~1.4 ms "
        "(FWHM 0.82 + tau600~0.53; NB joint-fit alpha railed at 6.0) -> ~4 el.",
    },
    "wilhelm": {
        "id": "253635173",
        "dm": 602.346,
        "fwhm_ms": 0.3887,
        "upchan": 64,
        "recoverable": True,
        "h5_relpath": "2022/12/03/astro_253635173/Run_UpdatedCalSep25/singlebeam_253635173.h5",
        "note": "MW-floor-dominated (25.9 kHz floor); U=64 ~4.25 ch (tightest freq margin of batch 2). "
        "Time: profile ~1.0 ms (tau600~1.04, alpha=2.71) -> ~3.2 el, marginal but passes.",
    },
    "oran": {
        "id": "224263996",
        "dm": 396.882,
        "fwhm_ms": 74.2,
        "upchan": 128,
        "recoverable": True,
        "h5_relpath": "2022/05/06/astro_224263996/Run_UpdatedCalSep25/singlebeam_224263996.h5",
        "note": "MW-floor-dominated (23.9 kHz floor); U=64 gives only 3.9 ch -> U=128 (3.05 kHz, "
        "~7.8 ch). 74 ms FWHM buys ~110 time el even at dt=0.655 ms. Mid-tier _upchannel path.",
    },
    "hamilton": {
        "id": "318353610",
        "dm": 518.799,
        "fwhm_ms": 0.2018,
        "upchan": 64,
        "recoverable": False,
        "h5_relpath": "2023/09/13/astro_318353610/singlebeam_318353610.h5",
        "note": "SINGLE-BLOCK MODE: floor 38.6 kHz needs U>=41 but profile ~0.3 ms (tau600~0.11; "
        "alpha railed at 6.0) spans <1 el at dt=0.328 ms. W*Dnu_d~12 so the spectrum is recoverable "
        "from the one on-burst block; check block-straddle from time0 before trusting the ACF.",
    },
    "johndoeII": {
        "id": "311723353",
        "dm": 696.506,
        "fwhm_ms": 0.6957,
        "upchan": 512,
        "recoverable": False,
        "h5_relpath": "2023/08/14/astro_311723353/singlebeam_311723353.h5",
        "note": "SINGLE-BLOCK MODE: floor 5.8 kHz needs U>=269 -> U=512 (0.763 kHz, ~7.6 ch) but "
        "dt=2.62 ms vs profile ~1.7 ms (tau600~1.7, alpha railed low 1.37) -> ~0.7 el. W*Dnu_d~10. "
        "SLOW python-loop path (mahi-class). Straddle check mandatory.",
    },
    "isha": {
        "id": "252069198",
        "dm": 411.568,
        "fwhm_ms": 1.8053,
        "upchan": 256,
        "recoverable": False,
        "h5_relpath": "2022/11/13/astro_252069198/singlebeam_252069198.h5",
        "note": "NOT cleanly resolvable: railed DSA input + dominant scale at/past the time-smearing wall. Upper-bound only.",
    },
}


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _validate_factor(value: int) -> int:
    factor = int(value)
    if factor not in SUPPORTED_UPCHANNEL_FACTORS:
        allowed = ", ".join(str(item) for item in SUPPORTED_UPCHANNEL_FACTORS)
        raise ValueError(f"unsupported upchannel factor {factor}; choose one of {allowed}")
    return factor


def _container_digest(identity: str) -> str:
    match = re.search(r"sha256:([0-9a-fA-F]{64})(?:$|[^0-9a-fA-F])", identity)
    if match is None:
        raise ValueError("container identity must include a pinned sha256:<64 hex> digest")
    return match.group(1).lower()


def _load_dm_provenance(
    path: str | Path,
    *,
    event: str,
    require_ratified: bool,
) -> dict:
    provenance_path = Path(path)
    try:
        payload = json.loads(provenance_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read DM provenance {provenance_path}: {exc}") from exc
    if payload.get("schema") != DM_PROVENANCE_SCHEMA:
        raise ValueError(f"unsupported DM provenance schema in {provenance_path}")
    if payload.get("event") != event or payload.get("instrument") != INSTRUMENT_ID:
        raise ValueError(f"DM provenance event/instrument mismatch for {event}")
    dm = payload.get("dm_pc_cm3")
    if not isinstance(dm, (int, float)) or not np.isfinite(dm) or dm <= 0:
        raise ValueError("DM provenance must contain a positive finite dm_pc_cm3")
    status = payload.get("ratification_status")
    if status not in {"candidate", "ratified"}:
        raise ValueError("DM ratification_status must be candidate or ratified")
    if require_ratified:
        missing = [
            field
            for field in ("ratified_by", "ratified_at_utc", "decision_record")
            if not payload.get(field)
        ]
        if status != "ratified" or missing:
            detail = f"; missing {', '.join(missing)}" if missing else ""
            raise ValueError(
                f"authoritative generation requires ratified DM provenance for {event}{detail}"
            )
        try:
            datetime.fromisoformat(str(payload["ratified_at_utc"]).replace("Z", "+00:00"))
        except ValueError as exc:
            raise ValueError("ratified_at_utc must be an ISO-8601 timestamp") from exc
    return {
        "path": str(provenance_path.resolve()),
        "sha256": _sha256(provenance_path),
        "schema": payload["schema"],
        "event": event,
        "instrument": INSTRUMENT_ID,
        "dm_pc_cm3": float(dm),
        "ratification_status": status,
        "ratified_by": payload.get("ratified_by"),
        "ratified_at_utc": payload.get("ratified_at_utc"),
        "decision_record": payload.get("decision_record"),
    }


def _inspect_h5_layout(path: str | Path) -> dict:
    import h5py  # noqa: PLC0415

    with h5py.File(path, "r") as handle:
        if "tiedbeam_baseband" not in handle:
            raise ValueError(f"{path}: missing tiedbeam_baseband dataset")
        shape = tuple(int(value) for value in handle["tiedbeam_baseband"].shape)
    if len(shape) != 3 or shape[1] != 2:
        raise ValueError(
            f"{path}: tiedbeam_baseband must have shape (coarse_channel, 2, native_time)"
        )
    return {
        "tiedbeam_baseband_shape": list(shape),
        "coarse_channels_present": shape[0],
        "polarizations": shape[1],
        "native_time_samples": shape[2],
    }


def _expected_geometry(
    *,
    upchannel_factor: int,
    source_layout: dict,
    fine_oversample: int | None,
) -> dict:
    factor = _validate_factor(upchannel_factor)
    oversample = 2 if fine_oversample is None else int(fine_oversample)
    fft_size = factor * oversample
    hop_samples = 2 * factor
    native_time = int(source_layout["native_time_samples"])
    if oversample == 2:
        time_bins = native_time // fft_size
    else:
        time_bins = 0 if native_time < fft_size else 1 + (native_time - fft_size) // hop_samples
    return {
        "df_mhz": CHIME_COARSE_DF_MHZ / factor,
        "dt_s": CHIME_NATIVE_DT_S * hop_samples,
        "shape": [
            CHIME_NOMINAL_COARSE_CHANNELS * factor,
            time_bins,
        ],
        "grid": {
            "nominal_coarse_channels": CHIME_NOMINAL_COARSE_CHANNELS,
            "nominal_fine_positions": CHIME_NOMINAL_COARSE_CHANNELS * factor,
            "expected_measured_fine_positions": (
                int(source_layout["coarse_channels_present"]) * factor
            ),
            "frequency_order": "ascending",
            "authoritative_coordinate": "integer fine_channel_id",
            "nominal_frequency_formula_mhz": (
                "800.1953125 - (fine_channel_id + 0.5) * "
                "(0.390625 / upchannel_factor)"
            ),
        },
        "channelizer": {
            "fft_size": fft_size,
            "oversample": oversample,
            "hop_samples": hop_samples,
        },
    }


def _software_identity() -> dict:
    worker = Path(__file__).resolve()
    try:
        commit = subprocess.run(
            ["git", "-C", str(worker.parent), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = None
    return {
        "worker_path": str(worker),
        "worker_sha256": _sha256(worker),
        "git_commit": commit,
        "version": "chime-upchannel-factor-ladder-v1",
    }


def _dm_identity(dm_pc_cm3: float, provenance_sha256: str) -> str:
    value = f"{dm_pc_cm3:.9f}".rstrip("0").rstrip(".").replace(".", "p")
    return f"dm-{value}-{provenance_sha256[:12]}"


def _product_directory(
    out_root: str | Path,
    *,
    event: str,
    dm_provenance: dict,
    upchannel_factor: int,
    software: dict,
    container_identity: str,
) -> Path:
    container = _container_digest(container_identity)
    return (
        Path(out_root)
        / "instrument-chime-frb"
        / f"event-{event}"
        / _dm_identity(dm_provenance["dm_pc_cm3"], dm_provenance["sha256"])
        / f"u{_validate_factor(upchannel_factor):04d}"
        / f"worker-{software['worker_sha256'][:12]}"
        / f"container-{container[:12]}"
    )


def _assert_product_path_available(path: str | Path) -> None:
    product_dir = Path(path)
    if product_dir.exists():
        raise FileExistsError(
            f"existing product collision; refusing to overwrite immutable path: {product_dir}"
        )
    parent = product_dir
    while not parent.exists() and parent != parent.parent:
        parent = parent.parent
    if not parent.is_dir():
        raise FileNotFoundError(f"no existing parent directory for output path {product_dir}")


def _planned_product_manifest(
    *,
    event: str,
    source_h5: str | Path,
    dm_provenance: dict,
    upchannel_factor: int,
    container_identity: str,
    out_root: str | Path,
    fine_window: str | None,
    fine_oversample: int | None,
    save_polarizations: bool,
    command: list[str],
) -> dict:
    factor = _validate_factor(upchannel_factor)
    source_path = Path(source_h5)
    if not source_path.is_file():
        raise FileNotFoundError(f"staged CHIME H5 input missing: {source_path}")
    software = _software_identity()
    source_layout = _inspect_h5_layout(source_path)
    expected = _expected_geometry(
        upchannel_factor=factor,
        source_layout=source_layout,
        fine_oversample=fine_oversample,
    )
    out_root = Path(out_root)
    if not out_root.is_dir() or not os.access(out_root, os.W_OK):
        raise ValueError(f"output root must be an existing writable directory: {out_root}")
    product_dir = _product_directory(
        out_root,
        event=event,
        dm_provenance=dm_provenance,
        upchannel_factor=factor,
        software=software,
        container_identity=container_identity,
    )
    _assert_product_path_available(product_dir)
    product_names = {
        "stokes_i": "stokes_i.npy",
        "acf_input": "acf_input.npz",
        "nominal_frequencies": "frequencies_mhz.npy",
        "package_frequencies": "package_frequencies_mhz.npy",
        "source_valid": "source_valid.npy",
    }
    if save_polarizations:
        product_names["polarizations"] = [
            "polarization_0.npy",
            "polarization_1.npy",
        ]
    ratified = dm_provenance["ratification_status"] == "ratified"
    return {
        "schema": PRODUCT_MANIFEST_SCHEMA,
        "status": "planned",
        "mode": "authoritative" if ratified else "preflight_only",
        "science_status": (
            "input_product_only_pending_owner_review"
            if ratified
            else "diagnostic_plan_blocked_pending_ratified_dm"
        ),
        "identity": {
            "event": event,
            "instrument": INSTRUMENT_ID,
            "dm_provenance": dm_provenance,
            "upchannel_factor": factor,
            "software": software,
            "container": {
                "identity": container_identity,
                "digest": f"sha256:{_container_digest(container_identity)}",
            },
        },
        "source": {
            "h5_path": str(source_path.resolve()),
            "h5_sha256": _sha256(source_path),
            **source_layout,
        },
        "expected": expected,
        "product_directory": str(product_dir),
        "products": {
            name: (
                [{"path": item, "sha256": None} for item in value]
                if isinstance(value, list)
                else {"path": value, "sha256": None}
            )
            for name, value in product_names.items()
        },
        "acf_contract": {
            "owner_approved_mask_identity": None,
            "full_to_compact_mapping_identity": None,
            "status": "required_before_acf",
            "required_lag_support_fields": [
                "first_fit_lag",
                "min_support_pairs",
                "min_support_fraction",
            ],
            "note": (
                "The downstream ACF provenance record must bind the exact source-valid, "
                "owner-approved mask, mapping, and product-manifest bytes."
            ),
        },
        "processing": {
            "command": shlex.join(command),
            "time_shift": None,
            "fine_window": fine_window or "rectangular",
            "fine_oversample": fine_oversample or 2,
        },
    }


def _local_h5_path(name: str, relpath: str) -> Path:
    """Return the canonical h17 path for one project's staged singlebeam file."""
    return Path(LOCAL_H5_DIR) / name / Path(relpath).name


def _detected_products(spec: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return Stokes-I and per-polarization power waterfalls.

    ``_upchannel`` returns ``(npol, ntime, nfreq)`` complex voltages.  Keep the
    two detected polarization streams separate until after validation: their
    receiver-noise realizations are the independent inputs required by the
    high-band cross-ACF experiment.  The historical Stokes-I product remains
    their exact sum.
    """
    values = np.asarray(spec)
    if values.ndim != 3 or values.shape[0] != 2:
        raise ValueError("upchannelized spectrum must have shape (2, ntime, nfreq)")
    per_pol = np.transpose(np.abs(values) ** 2, (0, 2, 1))
    return np.sum(per_pol, axis=0), per_pol


def _restore_nominal_fine_grid(
    stokes_i: np.ndarray,
    per_pol: np.ndarray,
    package_freq_mhz: np.ndarray,
    fine_channel_ids: np.ndarray,
    upchan_factor: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Scatter retained fine channels onto the nominal CHIME frequency grid.

    Integer channel identifiers are authoritative.  Missing measurements remain
    NaN and are also represented by an explicit Boolean validity mask.  Preserve
    the package-reported frequency coordinate separately: baseband-analysis 1.9.0
    uses an inclusive linspace whose values differ by at most half a fine channel
    from the exact nominal-bin centers.
    """
    values = np.asarray(stokes_i)
    polarizations = np.asarray(per_pol)
    package_freq = np.asarray(package_freq_mhz, dtype=np.float64)
    raw_ids = np.asarray(fine_channel_ids)
    if not isinstance(upchan_factor, (int, np.integer)) or upchan_factor < 1:
        raise ValueError("upchan_factor must be a positive integer")
    if raw_ids.ndim != 1 or not np.issubdtype(raw_ids.dtype, np.integer):
        raise ValueError("fine channel identifiers must be a 1-D integer array")
    ids = raw_ids.astype(np.int64, copy=False)
    if values.ndim != 2 or values.shape[0] != ids.size:
        raise ValueError("Stokes-I rows must match fine channel identifiers")
    if polarizations.ndim != 3 or polarizations.shape[1:] != values.shape:
        raise ValueError("polarization rows and time bins must match Stokes I")
    if package_freq.shape != ids.shape:
        raise ValueError("package frequency rows must match fine channel identifiers")
    if np.unique(ids).size != ids.size:
        raise ValueError("fine channel identifiers must be unique")

    total = CHIME_NOMINAL_COARSE_CHANNELS * int(upchan_factor)
    if ids.size and (ids.min() < 0 or ids.max() >= total):
        raise ValueError(f"fine channel identifiers outside nominal range 0:{total}")

    full_package_freq = np.linspace(
        CHIME_TOP_EDGE_MHZ, CHIME_BOTTOM_EDGE_MHZ, total
    )
    if not np.allclose(
        package_freq, full_package_freq[ids], rtol=0.0, atol=1e-10
    ):
        raise ValueError("package frequency coordinate does not match fine identifiers")
    fine_width_mhz = CHIME_COARSE_DF_MHZ / int(upchan_factor)
    nominal_freq = CHIME_TOP_EDGE_MHZ - (
        np.arange(total, dtype=np.float64) + 0.5
    ) * fine_width_mhz

    full_stokes = np.full((total, values.shape[1]), np.nan, dtype=values.dtype)
    full_per_pol = np.full(
        (polarizations.shape[0], total, values.shape[1]),
        np.nan,
        dtype=polarizations.dtype,
    )
    valid = np.zeros(total, dtype=bool)
    full_stokes[ids] = values
    full_per_pol[:, ids] = polarizations
    valid[ids] = True
    return full_stokes, full_per_pol, nominal_freq, full_package_freq, valid


def _waterfall(
    h5_path: str,
    dm: float,
    U: int,
    *,
    time_shift: bool = True,
    fine_window: str | None = None,
    fine_oversample: int | None = None,
):
    """Coherently-dedispersed, upchannelized Stokes-I waterfall (n_fine_freq, n_time) + freq[MHz].

    We assemble the chain by hand rather than via baseband_analysis.analysis.waterfall_from_beamformed
    because that function is BROKEN in this image (v1.9.0): it passes upchannel()'s 3-tuple return
    straight into incoherent_dedisp, which does `matrix_in.copy()` -> AttributeError on the tuple.
    The pieces we use ARE the package's: coherent_dedisp + the internal _upchannel. After coherent
    dedispersion the baseband is fully de-chirped, so no incoherent_dedisp step is needed -- we form
    Stokes I directly from the upchannelized complex spectrum.
    """
    from baseband_analysis.core.bbdata import BBData  # noqa: PLC0415
    from baseband_analysis.core.dedispersion import coherent_dedisp  # noqa: PLC0415
    from baseband_analysis.core.sampling import _upchannel  # noqa: PLC0415

    data = BBData.from_file(h5_path)
    time0 = data["time0"][:]
    source_metadata = {
        "delta_time": CHIME_NATIVE_DT_S,
        "fpga_count": np.asarray(time0["fpga_count"], dtype=np.uint64).tolist(),
        "freq_mhz": np.asarray(data.index_map["freq"]["centre"], dtype=float).tolist(),
        "freq_id": np.asarray(data.index_map["freq"]["id"], dtype=np.uint32).tolist(),
    }
    # coherent_dedisp returns the transformed array unless write=True.  Keep
    # that return value: upchannelizing data["tiedbeam_baseband"] would silently
    # use the original dispersed voltages.  For current products time_shift is
    # disabled because that operation is circular in short channel buffers;
    # padded alignment is applied downstream from time0 metadata.
    dedispersed = coherent_dedisp(data, dm, time_shift=time_shift)

    # _upchannel returns (spec, freq, chan_id): spec is (npol, nblock, nfine) complex, freq the
    # fine-channel centres (MHz) ordered high->low. upchan factor U = fftsize/downfreq.
    freq_id = data.index_map["freq"]["id"][:]
    if fine_window is None and fine_oversample is None:
        spec, freq, fine_channel_ids = _upchannel(
            dedispersed,
            freq_id=freq_id,
            fftsize=2 * U,
            downfreq=2,
        )
        channelizer_metadata = {
            "implementation": "baseband_analysis_1.9.0_private_upchannel",
            "window": "rectangular",
            "upchannel_factor": int(U),
            "oversample": 2,
            "fft_size": int(2 * U),
            "downfreq": 2,
            "hop_samples": int(2 * U),
            "frame_center_offset_samples": (2 * U - 1) / 2.0,
            "normalization": "package_exact",
        }
    elif fine_window is None or fine_oversample is None:
        raise ValueError("fine_window and fine_oversample must be supplied together")
    else:
        from windowed_upchan import windowed_upchannel  # noqa: PLC0415

        spec, freq, fine_channel_ids, channelizer_metadata = windowed_upchannel(
            dedispersed,
            freq_id,
            upchan_factor=U,
            window=fine_window,
            oversample=fine_oversample,
        )
    stokes_i, per_pol = _detected_products(spec)
    return (
        stokes_i,
        per_pol,
        np.asarray(freq, dtype=np.float64),
        np.asarray(fine_channel_ids, dtype=np.int64),
        source_metadata,
        channelizer_metadata,
    )


def recover_target(
    plan: dict,
    run_unresolvable: bool = False,
    save_polarizations: bool = False,
    time_shift: bool = True,
    fine_window: str | None = None,
    fine_oversample: int | None = None,
) -> Path:
    """Materialize one already-preflighted authoritative product."""
    name = plan["identity"]["event"]
    t = TARGETS[name]
    if not t["recoverable"] and not run_unresolvable:
        raise SystemExit(
            f"{name} is flagged NOT cleanly resolvable ({t['note']}). "
            f"Re-run with --run-unresolvable to produce a lower-confidence upper-bound spectrum."
        )
    dm_provenance = plan["identity"]["dm_provenance"]
    if (
        plan.get("mode") != "authoritative"
        or dm_provenance.get("ratification_status") != "ratified"
    ):
        raise ValueError("production accepts only an authoritative plan with a ratified DM")
    if _sha256(Path(__file__).resolve()) != plan["identity"]["software"]["worker_sha256"]:
        raise ValueError(f"{name}: worker bytes changed after preflight")
    h5_path = plan["source"]["h5_path"]
    if _sha256(h5_path) != plan["source"]["h5_sha256"]:
        raise ValueError(f"{name}: source H5 changed after preflight")
    U = _validate_factor(plan["identity"]["upchannel_factor"])
    product_dir = Path(plan["product_directory"])
    _assert_product_path_available(product_dir)
    (
        stokes_i,
        per_pol,
        package_freq,
        fine_channel_ids,
        source_metadata,
        channelizer_metadata,
    ) = _waterfall(
        h5_path,
        dm_provenance["dm_pc_cm3"],
        U,
        time_shift=time_shift,
        fine_window=fine_window,
        fine_oversample=fine_oversample,
    )

    (
        stokes_i,
        per_pol,
        freq,
        package_freq,
        source_valid,
    ) = _restore_nominal_fine_grid(
        stokes_i,
        per_pol,
        package_freq,
        fine_channel_ids,
        U,
    )

    # Ascending frequency to match the FLITS BurstDataset convention.
    if freq[0] > freq[-1]:
        freq = freq[::-1]
        package_freq = package_freq[::-1]
        stokes_i = stokes_i[::-1, :]
        per_pol = per_pol[:, ::-1, :]
        source_valid = source_valid[::-1]

    df_fine = CHIME_COARSE_DF_MHZ / U
    n_fine, n_time = stokes_i.shape

    # --- ponytail self-check: the recovered grid must match the requested upchannelization ---
    # NaN channels are EXPECTED (CHIME masks RFI/missing channels); the downstream ACF uses nansum.
    # So require a healthy finite FRACTION, not all-finite -- only an all-NaN/empty result is a failure.
    assert n_time > 0, f"{name}: empty time axis"
    expected_fine = CHIME_NOMINAL_COARSE_CHANNELS * U
    assert n_fine == expected_fine, (
        f"{name}: {n_fine} fine positions != nominal {expected_fine} (U={U})"
    )
    expected_measured = len(source_metadata["freq_id"]) * U
    assert int(source_valid.sum()) == expected_measured, (
        f"{name}: {source_valid.sum()} measured positions != expected {expected_measured}"
    )
    assert not np.isfinite(stokes_i[~source_valid]).any(), (
        f"{name}: missing source positions contain finite Stokes-I values"
    )
    finite_frac = float(np.isfinite(stokes_i).mean())
    assert finite_frac > 0.3, f"{name}: only {finite_frac:.1%} finite Stokes-I -- effectively empty"
    measured_df = abs(np.nanmedian(np.diff(freq)))
    assert np.isclose(measured_df, df_fine, rtol=0.05), (
        f"{name}: fine channel width {measured_df:.6f} MHz != expected {df_fine:.6f} MHz (U={U})"
    )
    expected_shape = tuple(plan["expected"]["shape"])
    if stokes_i.shape != expected_shape:
        raise ValueError(
            f"{name}: realized shape {stokes_i.shape} != preflight shape {expected_shape}"
        )

    product_dir.mkdir(parents=True, exist_ok=False)
    spec_path = product_dir / "stokes_i.npy"
    np.save(spec_path, stokes_i.astype(np.float32))
    polarization_paths = []
    if save_polarizations:
        for pol_index, power in enumerate(per_pol):
            path = product_dir / f"polarization_{pol_index}.npy"
            np.save(path, power.astype(np.float32))
            polarization_paths.append(path)
    frequency_path = product_dir / "frequencies_mhz.npy"
    np.save(frequency_path, freq)
    package_frequency_path = product_dir / "package_frequencies_mhz.npy"
    np.save(package_frequency_path, package_freq)
    source_valid_path = product_dir / "source_valid.npy"
    np.save(source_valid_path, source_valid)
    dt_s = CHIME_NATIVE_DT_S * int(channelizer_metadata["hop_samples"])
    times_s = np.arange(n_time, dtype=np.float64) * dt_s
    acf_input_path = product_dir / "acf_input.npz"
    np.savez(
        acf_input_path,
        power_2d=stokes_i[::-1, :].astype(np.float32),
        frequencies_mhz=freq,
        times_s=times_s,
    )

    metadata = json.loads(json.dumps(plan))
    metadata.update(
        {
            "status": "complete",
            "completed_at_utc": datetime.now(UTC).isoformat(),
            "source_time0": source_metadata,
            "channelizer": channelizer_metadata,
        }
    )
    metadata["processing"]["time_shift"] = bool(time_shift)
    metadata["realized"] = {
        "shape": [int(n_fine), int(n_time)],
        "df_mhz": float(measured_df),
        "dt_s": float(dt_s),
        "finite_fraction": finite_frac,
        "source_valid": {
            "measured_fine_positions": int(source_valid.sum()),
            "missing_fine_positions": int((~source_valid).sum()),
            "missing_value": "NaN",
        },
        "max_package_nominal_offset_mhz": float(
            np.max(np.abs(package_freq - freq))
        ),
    }
    metadata["products"] = {
        "stokes_i": {"path": spec_path.name, "sha256": _sha256(spec_path)},
        "acf_input": {
            "path": acf_input_path.name,
            "sha256": _sha256(acf_input_path),
        },
        "polarizations": [
            {"path": path.name, "sha256": _sha256(path)}
            for path in polarization_paths
        ],
        "nominal_frequencies": {
            "path": frequency_path.name,
            "sha256": _sha256(frequency_path),
        },
        "package_frequencies": {
            "path": package_frequency_path.name,
            "sha256": _sha256(package_frequency_path),
        },
        "source_valid": {
            "path": source_valid_path.name,
            "sha256": _sha256(source_valid_path),
        },
    }
    if fine_window is not None:
        companion = Path(__file__).with_name("windowed_upchan.py")
        metadata["channelizer_producer"] = {
            "path": str(companion.resolve()),
            "sha256": _sha256(companion),
        }
    metadata_path = product_dir / "manifest.json"
    with metadata_path.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(
        f"[{name}] U={U} shape={stokes_i.shape} df={df_fine * 1e3:.3f} kHz "
        f"dt={dt_s * 1e3:.4f} ms finite={finite_frac:.1%} -> {product_dir}"
    )
    return spec_path


def main(argv: list[str]) -> int:
    import argparse

    p = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    p.add_argument("targets", nargs="*", default=list(TARGETS), help="targets (default: all)")
    p.add_argument("--out", default=DEFAULT_OUT_DIR, help="immutable product-tree root")
    p.add_argument(
        "--upchannel-factor",
        type=int,
        choices=SUPPORTED_UPCHANNEL_FACTORS,
        action="append",
        help=(
            "factor to generate; repeat for a ladder. If omitted, use each target's "
            "historical candidate factor"
        ),
    )
    p.add_argument(
        "--dm-provenance-dir",
        type=Path,
        required=True,
        help="directory containing one <target>.json DM provenance artifact per target",
    )
    p.add_argument(
        "--container-identity",
        required=True,
        help="pinned container identity containing sha256:<64 hex>",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="preflight only; never upchannelize or write burst products",
    )
    p.add_argument(
        "--planned-manifest",
        type=Path,
        help="exclusive output path for the dry-run campaign plan",
    )
    p.add_argument(
        "--no-time-shift",
        action="store_true",
        help="skip coherent_dedisp time_shift (circular in short buffers); align downstream from time0 metadata",
    )
    p.add_argument(
        "--run-unresolvable",
        action="store_true",
        help="also process targets flagged NOT cleanly resolvable (isha, hamilton, johndoeII single-block), as an upper bound",
    )
    p.add_argument(
        "--save-polarizations",
        action="store_true",
        help="also retain separate detected polarization waterfalls for independent-noise tests",
    )
    p.add_argument(
        "--fine-window",
        choices=("rectangular", "hann", "blackmanharris"),
        help="use the local windowed channelizer (requires --fine-oversample)",
    )
    p.add_argument(
        "--fine-oversample",
        type=int,
        choices=(2, 4),
        help="FFT oversampling and complex-bin grouping (requires --fine-window)",
    )
    args = p.parse_args(argv)

    if (args.fine_window is None) != (args.fine_oversample is None):
        p.error("--fine-window and --fine-oversample must be supplied together")
    if args.dry_run != (args.planned_manifest is not None):
        p.error("--dry-run and --planned-manifest must be supplied together")

    targets = args.targets or list(TARGETS)
    unknown = [n for n in targets if n not in TARGETS]
    if unknown:
        raise SystemExit(f"unknown target(s) {unknown}; known: {list(TARGETS)}")
    blocked_targets = [
        name for name in targets if not TARGETS[name]["recoverable"]
    ]
    if blocked_targets and not args.dry_run and not args.run_unresolvable:
        raise SystemExit(
            "production preflight blocked for targets requiring --run-unresolvable: "
            + ", ".join(blocked_targets)
        )
    plans = []
    for name in targets:
        provenance = _load_dm_provenance(
            args.dm_provenance_dir / f"{name}.json",
            event=name,
            require_ratified=not args.dry_run,
        )
        factors = args.upchannel_factor or [TARGETS[name]["upchan"]]
        for factor in factors:
            source_h5 = _local_h5_path(name, TARGETS[name]["h5_relpath"])
            plan = _planned_product_manifest(
                event=name,
                source_h5=source_h5,
                dm_provenance=provenance,
                upchannel_factor=factor,
                container_identity=args.container_identity,
                out_root=args.out,
                fine_window=args.fine_window,
                fine_oversample=args.fine_oversample,
                save_polarizations=args.save_polarizations,
                command=[str(Path(__file__).resolve()), *argv],
            )
            plan["processing"]["time_shift"] = not args.no_time_shift
            plans.append(plan)

    product_directories = [item["product_directory"] for item in plans]
    if len(set(product_directories)) != len(product_directories):
        raise ValueError("duplicate target/factor product identity in requested campaign")

    if args.dry_run:
        plan_status = (
            "ready_for_authoritative_generation"
            if all(
                item["identity"]["dm_provenance"]["ratification_status"] == "ratified"
                for item in plans
            )
            else "blocked_pending_ratified_dm"
        )
        campaign = {
            "schema": PLAN_MANIFEST_SCHEMA,
            "status": plan_status,
            "science_status": "planning_only_no_burst_products_generated",
            "supported_upchannel_factors": list(SUPPORTED_UPCHANNEL_FACTORS),
            "planned_product_count": len(plans),
            "products": plans,
        }
        args.planned_manifest.parent.mkdir(parents=True, exist_ok=True)
        try:
            with args.planned_manifest.open("x", encoding="utf-8") as stream:
                stream.write(json.dumps(campaign, indent=2, sort_keys=True) + "\n")
        except FileExistsError as exc:
            raise FileExistsError(
                f"refusing to overwrite existing planned manifest: {args.planned_manifest}"
            ) from exc
        print(
            json.dumps(
                {
                    "status": plan_status,
                    "planned_product_count": len(plans),
                    "planned_manifest": str(args.planned_manifest),
                },
                sort_keys=True,
            )
        )
        return 0

    for plan in plans:
        recover_target(
            plan,
            run_unresolvable=args.run_unresolvable,
            save_polarizations=args.save_polarizations,
            time_shift=not args.no_time_shift,
            fine_window=args.fine_window,
            fine_oversample=args.fine_oversample,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
