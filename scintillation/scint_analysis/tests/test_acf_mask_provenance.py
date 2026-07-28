from __future__ import annotations

import json

import numpy as np
import pytest

from scint_analysis.acf_mask_provenance import (
    ProvenanceError,
    apply_verified_effective_mask,
    write_mapping_artifact,
)
from scint_analysis.core import DynamicSpectrum
from scint_analysis.freya_scintillation import prepare_spectrum_from_config
from scint_analysis.pipeline import ScintillationAnalysis


@pytest.mark.parametrize("upchannelization", [16, 64, 128, 256, 512])
@pytest.mark.parametrize("frequency_average", [1, 2, 4])
def test_full_grid_mapping_masks_before_frequency_averaging(
    tmp_path, upchannelization, frequency_average
):
    n_full = 1024 * upchannelization
    full_axis = np.arange(n_full, dtype=np.float64)
    source_valid = np.ones(n_full, dtype=bool)
    source_valid[frequency_average - 1 :: 97] = False
    compact_rows = np.flatnonzero(source_valid)
    compact_axis = full_axis[compact_rows]
    full_to_compact = np.full(n_full, -1, dtype=np.int64)
    full_to_compact[compact_rows] = np.arange(compact_rows.size)
    effective_mask = ~source_valid
    target_full = compact_rows[frequency_average]
    effective_mask[target_full] = True

    input_path = tmp_path / "burst.npz"
    np.savez(
        input_path,
        power_2d=np.ones((compact_rows.size, 2)),
        frequencies_mhz=compact_axis,
        times_s=np.arange(2),
    )
    source_valid_path = tmp_path / "source_valid.npy"
    effective_mask_path = tmp_path / "effective.npy"
    map_path = tmp_path / "owner-map.json"
    np.save(source_valid_path, source_valid)
    np.save(effective_mask_path, effective_mask)
    map_path.write_text('{"approval_status":"owner-approved"}\n')
    mapping_path = tmp_path / "mapping.npz"
    provenance_path = tmp_path / "mapping.json"
    write_mapping_artifact(
        mapping_path=mapping_path,
        provenance_path=provenance_path,
        input_path=input_path,
        full_frequency_axis=full_axis,
        compact_frequency_axis=compact_axis,
        source_valid_path=source_valid_path,
        map_path=map_path,
        effective_mask_path=effective_mask_path,
        full_to_compact=full_to_compact,
        event="zach",
        instrument="CHIME/FRB",
    )

    spectrum = DynamicSpectrum(
        np.ones((compact_rows.size, 2)), compact_axis, np.arange(2)
    )
    masked = apply_verified_effective_mask(
        spectrum,
        input_path=input_path,
        source_valid_path=source_valid_path,
        map_path=map_path,
        effective_mask_path=effective_mask_path,
        mapping_path=mapping_path,
        provenance_path=provenance_path,
        event="zach",
        instrument="CHIME/FRB",
    )
    compact_target = full_to_compact[target_full]
    assert np.all(np.ma.getmaskarray(masked.power)[compact_target])
    assert masked.power.count() == 2 * (compact_rows.size - 1)


def test_mapping_provenance_tamper_and_order_fail_closed(tmp_path):
    full_axis = np.arange(8, dtype=float)
    source_valid = np.array([1, 0, 1, 1, 0, 1, 1, 1], dtype=bool)
    rows = np.flatnonzero(source_valid)
    compact_axis = full_axis[rows]
    mapping = np.full(8, -1, dtype=np.int64)
    mapping[rows] = np.arange(rows.size)
    input_path = tmp_path / "burst.npz"
    np.savez(input_path, power_2d=np.ones((rows.size, 2)), frequencies_mhz=compact_axis, times_s=[0, 1])
    source_path = tmp_path / "valid.npy"
    mask_path = tmp_path / "mask.npy"
    map_path = tmp_path / "map.json"
    np.save(source_path, source_valid)
    np.save(mask_path, ~source_valid)
    map_path.write_text('{"approval_status":"owner-approved"}\n')
    mapping_path = tmp_path / "mapping.npz"
    provenance_path = tmp_path / "mapping.json"
    write_mapping_artifact(
        mapping_path=mapping_path,
        provenance_path=provenance_path,
        input_path=input_path,
        full_frequency_axis=full_axis,
        compact_frequency_axis=compact_axis,
        source_valid_path=source_path,
        map_path=map_path,
        effective_mask_path=mask_path,
        full_to_compact=mapping,
        event="zach",
        instrument="CHIME/FRB",
    )
    record = json.loads(provenance_path.read_text())
    record["event"] = "wrong"
    provenance_path.write_text(json.dumps(record))
    spectrum = DynamicSpectrum(np.ones((rows.size, 2)), compact_axis, np.arange(2))
    with pytest.raises(ProvenanceError, match="event"):
        apply_verified_effective_mask(
            spectrum,
            input_path=input_path,
            source_valid_path=source_path,
            map_path=map_path,
            effective_mask_path=mask_path,
            mapping_path=mapping_path,
            provenance_path=provenance_path,
            event="zach",
            instrument="CHIME/FRB",
        )


def _configured_case(tmp_path):
    full_axis = np.arange(32, dtype=float)
    source_valid = np.ones(32, dtype=bool)
    compact_axis = full_axis.copy()
    mapping = np.arange(32, dtype=np.int64)
    effective = np.zeros(32, dtype=bool)
    effective[[4, 8, 9]] = True
    power = np.ones((32, 4))
    # from_numpy_file reverses the stored power rows; this lands on compact row 4.
    power[-5] = 1e300
    input_path = tmp_path / "input.npz"
    source_path = tmp_path / "source.npy"
    effective_path = tmp_path / "effective.npy"
    map_path = tmp_path / "map.json"
    mapping_path = tmp_path / "mapping.npz"
    provenance_path = tmp_path / "provenance.json"
    np.savez(input_path, power_2d=power, frequencies_mhz=compact_axis, times_s=np.arange(4))
    np.save(source_path, source_valid)
    np.save(effective_path, effective)
    map_path.write_text('{"approval_status":"owner-approved"}\n')
    record = write_mapping_artifact(
        mapping_path=mapping_path,
        provenance_path=provenance_path,
        input_path=input_path,
        full_frequency_axis=full_axis,
        compact_frequency_axis=compact_axis,
        source_valid_path=source_path,
        map_path=map_path,
        effective_mask_path=effective_path,
        full_to_compact=mapping,
        event="zach",
        instrument="CHIME/FRB",
    )
    config = {
        "burst_id": "zach",
        "input_data_path": str(input_path),
        "analysis": {
            "bad_channel_mask": {
                "required": True,
                "source_valid_path": str(source_path),
                "map_path": str(map_path),
                "effective_mask_path": str(effective_path),
                "mapping_path": str(mapping_path),
                "provenance_path": str(provenance_path),
                "event": "zach",
                "instrument": "CHIME/FRB",
                "expected_hashes": {
                    name: record["sha256"][name]
                    for name in ("input", "source_valid", "owner_map", "effective_mask", "mapping")
                },
            },
            "grid_regularization": {"enable": False},
            "rfi_masking": {"disable": True},
        },
        "pipeline_options": {
            "downsample": {"f_factor": 2, "t_factor": 1},
            "save_intermediate_steps": False,
        },
    }
    return config, provenance_path


@pytest.mark.parametrize("entrypoint", ["function", "object"])
def test_both_preparation_entrypoints_apply_mask_before_regularization(
    monkeypatch, tmp_path, entrypoint
):
    config, _ = _configured_case(tmp_path)
    import scint_analysis.freya_scintillation as freya

    observed = []

    def inspect_before_regularization(spectrum, _config):
        row_mask = np.all(np.ma.getmaskarray(spectrum.power), axis=1)
        observed.append(row_mask.copy())
        assert row_mask[4] and row_mask[8] and row_mask[9]
        assert not np.any(np.isfinite(spectrum.power.compressed()) & (spectrum.power.compressed() > 1e200))
        return spectrum

    monkeypatch.setattr(freya, "apply_grid_regularization", inspect_before_regularization)
    monkeypatch.setattr(DynamicSpectrum, "mask_rfi", lambda self, _cfg: self)
    if entrypoint == "function":
        prepare_spectrum_from_config(config)
    else:
        ScintillationAnalysis(config).prepare_data()
    assert len(observed) == 1


def test_partial_and_all_masked_frequency_blocks(monkeypatch, tmp_path):
    config, _ = _configured_case(tmp_path)
    monkeypatch.setattr(DynamicSpectrum, "mask_rfi", lambda self, _cfg: self)
    spectrum, *_ = prepare_spectrum_from_config(config)
    row_mask = np.all(np.ma.getmaskarray(spectrum.power), axis=1)
    # factor-2 block [4,5] retains row 5; block [8,9] is fully masked.
    assert not row_mask[2]
    assert row_mask[4]
    assert np.allclose(spectrum.power[2].compressed(), 1.0)


def test_cache_identity_tracks_provenance_bytes(tmp_path):
    config, provenance_path = _configured_case(tmp_path)
    first = ScintillationAnalysis(config)._get_cache_path("processed_spectrum")
    provenance_path.write_text(provenance_path.read_text() + "\n")
    second = ScintillationAnalysis(config)._get_cache_path("processed_spectrum")
    assert first != second


def test_required_mask_missing_provenance_fails_closed(tmp_path):
    config, provenance_path = _configured_case(tmp_path)
    provenance_path.unlink()
    spectrum = DynamicSpectrum(np.ones((32, 2)), np.arange(32.0), np.arange(2))
    from scint_analysis.acf_mask_provenance import apply_configured_effective_mask

    with pytest.raises(ProvenanceError):
        apply_configured_effective_mask(spectrum, config)


@pytest.mark.parametrize("failure", ["axis-order", "duplicate-row"])
def test_mapping_writer_rejects_axis_or_row_order_failure(tmp_path, failure):
    full_axis = np.arange(6, dtype=float)
    compact_axis = full_axis.copy()
    mapping = np.arange(6, dtype=np.int64)
    if failure == "axis-order":
        compact_axis[[1, 2]] = compact_axis[[2, 1]]
    else:
        mapping[2] = 1
    input_path = tmp_path / "input.npz"
    source_path = tmp_path / "source.npy"
    effective_path = tmp_path / "effective.npy"
    map_path = tmp_path / "map.json"
    np.savez(input_path, power_2d=np.ones((6, 2)), frequencies_mhz=compact_axis, times_s=[0, 1])
    np.save(source_path, np.ones(6, dtype=bool))
    np.save(effective_path, np.zeros(6, dtype=bool))
    map_path.write_text("{}")
    with pytest.raises(ProvenanceError):
        write_mapping_artifact(
            mapping_path=tmp_path / "mapping.npz",
            provenance_path=tmp_path / "provenance.json",
            input_path=input_path,
            full_frequency_axis=full_axis,
            compact_frequency_axis=compact_axis,
            source_valid_path=source_path,
            map_path=map_path,
            effective_mask_path=effective_path,
            full_to_compact=mapping,
            event="zach",
            instrument="CHIME/FRB",
        )
