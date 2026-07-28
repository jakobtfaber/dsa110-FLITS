"""Unit checks for provenance-preserving CHIME detected products."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.fft import fft, fftshift

MODULE = Path(__file__).with_name("upchannelize_chime.py")
WINDOWED_MODULE = Path(__file__).with_name("windowed_upchan.py")


def _module():
    spec = importlib.util.spec_from_file_location("upchannelize_chime", MODULE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _windowed_module():
    spec = importlib.util.spec_from_file_location("windowed_upchan", WINDOWED_MODULE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _package_reference(wfall, freq_id, *, fftsize, downfreq):
    """Dependency-light copy of baseband_analysis 1.9.0 ``_upchannel``."""
    values = np.swapaxes(np.swapaxes(np.asarray(wfall), 0, 1), 1, 2)
    npol, nsamp, nchan = values.shape
    upchan = fftsize // downfreq
    nblock = nsamp // fftsize
    spectrum = np.zeros((npol, nblock, nchan * upchan), dtype=np.complex64)
    channel_ids = np.zeros(nchan * upchan, dtype=int)
    full_band = np.linspace(800.1953125, 400.1953125, upchan * 1024)
    for pol in range(npol):
        for block in range(nblock):
            for channel in range(nchan):
                time_series = values[
                    pol, block * fftsize : (block + 1) * fftsize, channel
                ].copy()
                transformed = fftshift(fft(time_series))
                transformed = transformed.reshape(upchan, downfreq).mean(axis=1).copy()
                start = channel * upchan
                spectrum[pol, block, start : start + upchan] = transformed
                channel_ids[start : start + upchan] = np.arange(
                    upchan * freq_id[channel], upchan * freq_id[channel] + upchan
                )
    return spectrum, full_band[channel_ids], channel_ids


def test_detected_products_preserve_independent_polarizations_and_stokes_sum():
    module = _module()
    rng = np.random.default_rng(20260714)
    voltages = rng.normal(size=(2, 7, 11)) + 1j * rng.normal(size=(2, 7, 11))

    stokes_i, per_pol = module._detected_products(voltages)

    assert per_pol.shape == (2, 11, 7)
    np.testing.assert_allclose(per_pol[0], np.abs(voltages[0]).T ** 2)
    np.testing.assert_allclose(per_pol[1], np.abs(voltages[1]).T ** 2)
    np.testing.assert_allclose(stokes_i, per_pol.sum(axis=0))


def test_detected_products_rejects_missing_polarization_axis():
    module = _module()

    with np.testing.assert_raises_regex(ValueError, "shape"):
        module._detected_products(np.ones((8, 16), dtype=complex))


def test_nominal_grid_restoration_preserves_measured_values_and_masks_gaps():
    module = _module()
    upchan = 2
    fine_ids = np.array([2, 3, 8, 9])
    package_grid = np.linspace(800.1953125, 400.1953125, 1024 * upchan)
    package_freq = package_grid[fine_ids]
    stokes = np.arange(12, dtype=np.float32).reshape(4, 3)
    per_pol = np.stack((stokes, stokes + 100.0))

    restored = module._restore_nominal_fine_grid(
        stokes, per_pol, package_freq, fine_ids, upchan
    )
    full_stokes, full_per_pol, nominal_freq, full_package_freq, valid = restored

    assert full_stokes.shape == (2048, 3)
    assert full_per_pol.shape == (2, 2048, 3)
    assert valid.shape == (2048,)
    assert valid.sum() == 4
    np.testing.assert_array_equal(full_stokes[fine_ids], stokes)
    np.testing.assert_array_equal(full_per_pol[:, fine_ids], per_pol)
    assert np.isnan(full_stokes[~valid]).all()
    assert np.isnan(full_per_pol[:, ~valid]).all()
    np.testing.assert_array_equal(full_package_freq, package_grid)
    expected_nominal = 800.1953125 - (np.arange(2048) + 0.5) * (0.390625 / upchan)
    np.testing.assert_array_equal(nominal_freq, expected_nominal)


def test_nominal_grid_restoration_rejects_bad_fine_identifiers():
    module = _module()
    stokes = np.ones((2, 3))
    per_pol = np.ones((2, 2, 3))
    package_grid = np.linspace(800.1953125, 400.1953125, 2048)

    with np.testing.assert_raises_regex(ValueError, "unique"):
        module._restore_nominal_fine_grid(
            stokes, per_pol, package_grid[[2, 2]], np.array([2, 2]), 2
        )
    with np.testing.assert_raises_regex(ValueError, "range"):
        module._restore_nominal_fine_grid(
            stokes, per_pol, np.array([1.0, 2.0]), np.array([2, 2048]), 2
        )


def test_nominal_grid_restoration_checks_package_frequency_mapping():
    module = _module()
    stokes = np.ones((2, 3))
    per_pol = np.ones((2, 2, 3))

    with np.testing.assert_raises_regex(ValueError, "package frequency"):
        module._restore_nominal_fine_grid(
            stokes, per_pol, np.array([700.0, 600.0]), np.array([2, 3]), 2
        )


def _dm_artifact(path: Path, *, status: str = "candidate") -> Path:
    payload = {
        "schema": "faber2026-chime-dm-provenance-v1",
        "event": "zach",
        "instrument": "CHIME/FRB",
        "dm_pc_cm3": 262.368,
        "ratification_status": status,
    }
    if status == "ratified":
        payload.update(
            {
                "ratified_by": "owner",
                "ratified_at_utc": "2026-07-29T00:00:00Z",
                "decision_record": "owner-queue:dm-zach",
            }
        )
    path.write_text(json.dumps(payload))
    return path


def test_factor_ladder_includes_u32_and_rejects_other_values():
    module = _module()

    assert module.SUPPORTED_UPCHANNEL_FACTORS == (16, 32, 64, 128, 256, 512)
    assert module._validate_factor(32) == 32
    with pytest.raises(ValueError, match="unsupported upchannel factor"):
        module._validate_factor(17)


@pytest.mark.parametrize("factor", [16, 32, 64, 128, 256, 512])
def test_factor_specific_resolution_and_shape(factor):
    module = _module()
    geometry = module._expected_geometry(
        upchannel_factor=factor,
        source_layout={
            "coarse_channels_present": 768,
            "native_time_samples": 4096,
        },
        fine_oversample=None,
    )

    assert geometry["df_mhz"] == module.CHIME_COARSE_DF_MHZ / factor
    assert geometry["dt_s"] == module.CHIME_NATIVE_DT_S * 2 * factor
    assert geometry["shape"] == [1024 * factor, 4096 // (2 * factor)]
    assert geometry["grid"]["expected_measured_fine_positions"] == 768 * factor


def test_authoritative_dm_gate_requires_ratification_receipt(tmp_path):
    module = _module()
    candidate = _dm_artifact(tmp_path / "zach.json")

    with pytest.raises(ValueError, match="requires ratified DM provenance"):
        module._load_dm_provenance(
            candidate,
            event="zach",
            require_ratified=True,
        )

    ratified = _dm_artifact(tmp_path / "zach-ratified.json", status="ratified")
    loaded = module._load_dm_provenance(
        ratified,
        event="zach",
        require_ratified=True,
    )
    assert loaded["ratification_status"] == "ratified"
    assert loaded["sha256"]


def test_immutable_product_identity_refuses_collision(tmp_path):
    module = _module()
    software = {"worker_sha256": "1" * 64}
    dm = {"dm_pc_cm3": 262.368, "sha256": "2" * 64}
    product_dir = module._product_directory(
        tmp_path,
        event="zach",
        dm_provenance=dm,
        upchannel_factor=32,
        software=software,
        container_identity=f"image@sha256:{'3' * 64}",
    )
    assert "event-zach" in str(product_dir)
    assert "u0032" in str(product_dir)
    assert "worker-111111111111" in str(product_dir)
    product_dir.mkdir(parents=True)
    with pytest.raises(FileExistsError, match="refusing to overwrite immutable path"):
        module._assert_product_path_available(product_dir)


def test_planned_manifest_contract_is_complete(tmp_path, monkeypatch):
    module = _module()
    source = tmp_path / "singlebeam.h5"
    source.write_bytes(b"fixture")
    dm_path = _dm_artifact(tmp_path / "zach.json")
    dm = module._load_dm_provenance(
        dm_path,
        event="zach",
        require_ratified=False,
    )
    monkeypatch.setattr(
        module,
        "_software_identity",
        lambda: {
            "worker_path": "/work/upchannelize_chime.py",
            "worker_sha256": "1" * 64,
            "git_commit": "a" * 40,
            "version": "test",
        },
    )
    monkeypatch.setattr(
        module,
        "_inspect_h5_layout",
        lambda _path: {
            "tiedbeam_baseband_shape": [768, 2, 4096],
            "coarse_channels_present": 768,
            "polarizations": 2,
            "native_time_samples": 4096,
        },
    )
    plan = module._planned_product_manifest(
        event="zach",
        source_h5=source,
        dm_provenance=dm,
        upchannel_factor=32,
        container_identity=f"image@sha256:{'3' * 64}",
        out_root=tmp_path,
        fine_window=None,
        fine_oversample=None,
        save_polarizations=True,
        command=["python", "upchannelize_chime.py", "zach"],
    )

    assert plan["schema"] == module.PRODUCT_MANIFEST_SCHEMA
    assert plan["science_status"] == "diagnostic_plan_blocked_pending_ratified_dm"
    assert plan["identity"]["upchannel_factor"] == 32
    assert plan["source"]["h5_sha256"]
    assert plan["expected"]["df_mhz"] == module.CHIME_COARSE_DF_MHZ / 32
    assert plan["expected"]["dt_s"] == module.CHIME_NATIVE_DT_S * 64
    assert plan["expected"]["shape"] == [32768, 64]
    assert plan["expected"]["grid"]["expected_measured_fine_positions"] == 24576
    assert plan["products"]["source_valid"]["path"] == "source_valid.npy"
    assert len(plan["products"]["polarizations"]) == 2
    assert plan["acf_contract"]["status"] == "required_before_acf"
    assert set(plan["acf_contract"]["required_lag_support_fields"]) == {
        "first_fit_lag",
        "min_support_pairs",
        "min_support_fraction",
    }
    assert plan["processing"]["command"]


def test_dry_run_emits_plan_without_generating_product(tmp_path, monkeypatch):
    module = _module()
    h5_root = tmp_path / "h5"
    source = h5_root / "zach" / "singlebeam_210456524.h5"
    source.parent.mkdir(parents=True)
    source.write_bytes(b"fixture")
    dm_dir = tmp_path / "dm"
    dm_dir.mkdir()
    _dm_artifact(dm_dir / "zach.json")
    out_root = tmp_path / "products"
    out_root.mkdir()
    plan_path = tmp_path / "planned.json"
    monkeypatch.setattr(module, "LOCAL_H5_DIR", str(h5_root))
    monkeypatch.setattr(
        module,
        "_software_identity",
        lambda: {
            "worker_path": "/work/upchannelize_chime.py",
            "worker_sha256": "1" * 64,
            "git_commit": "a" * 40,
            "version": "test",
        },
    )
    monkeypatch.setattr(
        module,
        "_inspect_h5_layout",
        lambda _path: {
            "tiedbeam_baseband_shape": [768, 2, 4096],
            "coarse_channels_present": 768,
            "polarizations": 2,
            "native_time_samples": 4096,
        },
    )
    result = module.main(
        [
            "zach",
            "--upchannel-factor",
            "32",
            "--dm-provenance-dir",
            str(dm_dir),
            "--container-identity",
            f"image@sha256:{'3' * 64}",
            "--out",
            str(out_root),
            "--dry-run",
            "--planned-manifest",
            str(plan_path),
        ]
    )

    assert result == 0
    campaign = json.loads(plan_path.read_text())
    assert campaign["status"] == "blocked_pending_ratified_dm"
    assert campaign["planned_product_count"] == 1
    product_dir = Path(campaign["products"][0]["product_directory"])
    assert not product_dir.exists()
    assert list(out_root.iterdir()) == []


def test_completed_receipt_hashes_products_and_preserves_expected_geometry(
    tmp_path, monkeypatch
):
    module = _module()
    source = tmp_path / "singlebeam.h5"
    source.write_bytes(b"fixture")
    dm_path = _dm_artifact(tmp_path / "zach.json", status="ratified")
    dm = module._load_dm_provenance(
        dm_path,
        event="zach",
        require_ratified=True,
    )
    worker_sha = module._sha256(module.Path(module.__file__).resolve())
    monkeypatch.setattr(
        module,
        "_software_identity",
        lambda: {
            "worker_path": module.__file__,
            "worker_sha256": worker_sha,
            "git_commit": "a" * 40,
            "version": "test",
        },
    )
    monkeypatch.setattr(
        module,
        "_inspect_h5_layout",
        lambda _path: {
            "tiedbeam_baseband_shape": [512, 2, 64],
            "coarse_channels_present": 512,
            "polarizations": 2,
            "native_time_samples": 64,
        },
    )
    out_root = tmp_path / "products"
    out_root.mkdir()
    plan = module._planned_product_manifest(
        event="zach",
        source_h5=source,
        dm_provenance=dm,
        upchannel_factor=16,
        container_identity=f"image@sha256:{'3' * 64}",
        out_root=out_root,
        fine_window=None,
        fine_oversample=None,
        save_polarizations=False,
        command=["python", "upchannelize_chime.py", "zach"],
    )
    fine_ids = np.arange(512 * 16, dtype=np.int64)
    package_grid = np.linspace(800.1953125, 400.1953125, 1024 * 16)
    stokes = np.ones((fine_ids.size, 2), dtype=np.float32)
    per_pol = np.stack((stokes, stokes))
    monkeypatch.setattr(
        module,
        "_waterfall",
        lambda *_args, **_kwargs: (
            stokes,
            per_pol,
            package_grid[fine_ids],
            fine_ids,
            {
                "delta_time": module.CHIME_NATIVE_DT_S,
                "fpga_count": [0] * 512,
                "freq_mhz": [600.0] * 512,
                "freq_id": list(range(512)),
            },
            {
                "implementation": "test",
                "window": "rectangular",
                "upchannel_factor": 16,
                "oversample": 2,
                "fft_size": 32,
                "downfreq": 2,
                "hop_samples": 32,
            },
        ),
    )

    stokes_path = module.recover_target(plan, time_shift=False)
    manifest_path = stokes_path.parent / "manifest.json"
    manifest = json.loads(manifest_path.read_text())

    assert manifest["status"] == "complete"
    assert manifest["mode"] == "authoritative"
    assert manifest["identity"]["dm_provenance"]["ratification_status"] == "ratified"
    assert manifest["expected"]["shape"] == [16384, 2]
    assert manifest["realized"]["shape"] == [16384, 2]
    for name in (
        "stokes_i",
        "acf_input",
        "nominal_frequencies",
        "package_frequencies",
        "source_valid",
    ):
        entry = manifest["products"][name]
        assert entry["sha256"] == module._sha256(stokes_path.parent / entry["path"])
    assert manifest["acf_contract"]["owner_approved_mask_identity"] is None
    assert manifest["acf_contract"]["full_to_compact_mapping_identity"] is None


def test_rectangular_oversample_two_is_package_equivalent():
    module = _windowed_module()
    rng = np.random.default_rng(20260714)
    wfall = (
        rng.normal(size=(3, 2, 47)) + 1j * rng.normal(size=(3, 2, 47))
    ).astype(np.complex64)
    freq_id = np.array([4, 8, 11])

    actual = module.windowed_upchannel(
        wfall,
        freq_id,
        upchan_factor=4,
        window="rectangular",
        oversample=2,
    )
    expected = _package_reference(wfall, freq_id, fftsize=8, downfreq=2)

    np.testing.assert_array_equal(actual[0], expected[0])
    np.testing.assert_array_equal(actual[1], expected[1])
    np.testing.assert_array_equal(actual[2], expected[2])
    assert actual[3]["hop_samples"] == 8
    assert actual[3]["frame_center_offset_samples"] == 3.5
    assert actual[3]["normalization"] == "package_exact"


def test_oversample_four_preserves_output_cadence_and_frequency_grid():
    module = _windowed_module()
    rng = np.random.default_rng(18)
    wfall = rng.normal(size=(2, 2, 64)) + 1j * rng.normal(size=(2, 2, 64))
    freq_id = np.array([2, 9])

    two = module.windowed_upchannel(
        wfall, freq_id, upchan_factor=4, window="hann", oversample=2
    )
    four = module.windowed_upchannel(
        wfall, freq_id, upchan_factor=4, window="hann", oversample=4
    )

    assert two[0].shape == (2, 8, 8)
    assert four[0].shape == (2, 7, 8)
    np.testing.assert_array_equal(four[1], two[1])
    np.testing.assert_array_equal(four[2], two[2])
    assert four[3]["fft_size"] == 16
    assert four[3]["downfreq"] == 4
    assert four[3]["hop_samples"] == 8
    assert four[3]["frame_center_offset_samples"] == 7.5


def test_exact_grouped_bin_normalization_preserves_white_noise_power():
    module = _windowed_module()
    rng = np.random.default_rng(1024)
    wfall = (
        rng.normal(size=(1, 2, 131072)) + 1j * rng.normal(size=(1, 2, 131072))
    ).astype(np.complex64)
    freq_id = np.array([5])
    powers = {}
    for window, oversample in (
        ("rectangular", 2),
        ("hann", 2),
        ("hann", 4),
        ("blackmanharris", 2),
        ("blackmanharris", 4),
    ):
        spectrum, _, _, metadata = module.windowed_upchannel(
            wfall,
            freq_id,
            upchan_factor=16,
            window=window,
            oversample=oversample,
        )
        powers[(window, oversample)] = float(np.mean(np.abs(spectrum) ** 2))
        assert np.isclose(metadata["grouped_noise_gain"], 16.0, rtol=1e-12)
    baseline = powers[("rectangular", 2)]
    for value in powers.values():
        assert np.isclose(value, baseline, rtol=0.03)


def test_windows_reduce_fractional_bin_far_sidelobes():
    module = _windowed_module()
    upchan = 64
    size = 2 * upchan
    samples = np.arange(size)
    tone = np.exp(2j * np.pi * 9.37 * samples / size)[None, None, :]
    freq_id = np.array([0])

    spectra = {}
    for window in ("rectangular", "hann", "blackmanharris"):
        spectrum, _, _, _ = module.windowed_upchannel(
            tone,
            freq_id,
            upchan_factor=upchan,
            window=window,
            oversample=2,
        )
        power = np.abs(spectrum[0, 0]) ** 2
        power /= power.max()
        spectra[window] = power

    peak = int(np.argmax(spectra["rectangular"]))
    far = np.ones(upchan, dtype=bool)
    far[max(0, peak - 3) : min(upchan, peak + 4)] = False
    assert spectra["hann"][far].max() < spectra["rectangular"][far].max() / 5
    assert spectra["blackmanharris"][far].max() < spectra["rectangular"][far].max() / 20
