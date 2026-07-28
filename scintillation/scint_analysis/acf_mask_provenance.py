"""Fail-closed provenance for full-grid RFI masks applied to compact spectra."""

from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path

import numpy as np

from .core import DynamicSpectrum


class ProvenanceError(ValueError):
    """The configured mask cannot be proven to belong to this spectrum."""


def _file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _array_sha256(value: np.ndarray) -> str:
    buffer = io.BytesIO()
    np.save(buffer, np.asarray(value), allow_pickle=False)
    return hashlib.sha256(buffer.getvalue()).hexdigest()


def write_mapping_artifact(
    *,
    mapping_path: str | Path,
    provenance_path: str | Path,
    input_path: str | Path,
    full_frequency_axis: np.ndarray,
    compact_frequency_axis: np.ndarray,
    source_valid_path: str | Path,
    map_path: str | Path,
    effective_mask_path: str | Path,
    full_to_compact: np.ndarray,
    event: str,
    instrument: str,
) -> dict:
    """Write an index mapping and a byte-bound provenance record.

    Mapping is explicit integer indexing. Frequency values are validation
    evidence only; they are never used to infer row correspondence.
    """
    full_axis = np.asarray(full_frequency_axis, dtype=np.float64)
    compact_axis = np.asarray(compact_frequency_axis, dtype=np.float64)
    mapping = np.asarray(full_to_compact, dtype=np.int64)
    source_valid = np.asarray(np.load(source_valid_path, allow_pickle=False), dtype=bool)
    effective_mask = np.asarray(np.load(effective_mask_path, allow_pickle=False), dtype=bool)
    _validate_mapping_arrays(
        full_axis, compact_axis, source_valid, effective_mask, mapping
    )
    mapping_path = Path(mapping_path)
    np.savez_compressed(
        mapping_path,
        full_frequency_mhz=full_axis,
        compact_frequency_mhz=compact_axis,
        full_to_compact=mapping,
    )
    record = {
        "schema": "faber2026-acf-full-to-compact-v1",
        "event": str(event),
        "instrument": str(instrument),
        "ordering": "full-grid-index-to-compact-row",
        "rows": {"full": int(full_axis.size), "compact": int(compact_axis.size)},
        "sha256": {
            "input": _file_sha256(input_path),
            "source_valid": _file_sha256(source_valid_path),
            "owner_map": _file_sha256(map_path),
            "effective_mask": _file_sha256(effective_mask_path),
            "mapping": _file_sha256(mapping_path),
            "full_frequency_axis": _array_sha256(full_axis),
            "compact_frequency_axis": _array_sha256(compact_axis),
        },
    }
    Path(provenance_path).write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return record


def _validate_mapping_arrays(
    full_axis: np.ndarray,
    compact_axis: np.ndarray,
    source_valid: np.ndarray,
    effective_mask: np.ndarray,
    mapping: np.ndarray,
) -> None:
    n_full = full_axis.size
    if any(value.ndim != 1 for value in (full_axis, compact_axis, source_valid, effective_mask, mapping)):
        raise ProvenanceError("mapping inputs must be one-dimensional")
    if not (source_valid.size == effective_mask.size == mapping.size == n_full):
        raise ProvenanceError("full-grid row counts disagree")
    if not np.all(np.isfinite(full_axis)) or not np.all(np.isfinite(compact_axis)):
        raise ProvenanceError("frequency axis contains non-finite values")
    mapped = mapping >= 0
    if not np.array_equal(mapped, source_valid):
        raise ProvenanceError("mapping presence does not exactly match source-valid rows")
    compact_ids = mapping[mapped]
    if not np.array_equal(np.sort(compact_ids), np.arange(compact_axis.size)):
        raise ProvenanceError("compact rows are missing, duplicated, or out of range")
    full_rows_by_compact = np.empty(compact_axis.size, dtype=np.int64)
    full_rows_by_compact[compact_ids] = np.flatnonzero(mapped)
    if not np.array_equal(compact_axis, full_axis[full_rows_by_compact]):
        raise ProvenanceError("frequency axis/order does not match explicit mapping")


def apply_verified_effective_mask(
    spectrum: DynamicSpectrum,
    *,
    input_path: str | Path,
    source_valid_path: str | Path,
    map_path: str | Path,
    effective_mask_path: str | Path,
    mapping_path: str | Path,
    provenance_path: str | Path,
    event: str,
    instrument: str,
    expected_hashes: dict[str, str] | None = None,
) -> DynamicSpectrum:
    """Validate every identity and mask compact rows in place."""
    try:
        record = json.loads(Path(provenance_path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ProvenanceError(f"cannot read mapping provenance: {exc}") from exc
    if record.get("schema") != "faber2026-acf-full-to-compact-v1":
        raise ProvenanceError("unsupported mapping provenance schema")
    for field, expected in (("event", event), ("instrument", instrument)):
        if record.get(field) != expected:
            raise ProvenanceError(f"{field} mismatch")
    if record.get("ordering") != "full-grid-index-to-compact-row":
        raise ProvenanceError("mapping ordering mismatch")

    paths = {
        "input": input_path,
        "source_valid": source_valid_path,
        "owner_map": map_path,
        "effective_mask": effective_mask_path,
        "mapping": mapping_path,
    }
    recorded_hashes = record.get("sha256", {})
    for name, path in paths.items():
        actual = _file_sha256(path)
        if recorded_hashes.get(name) != actual:
            raise ProvenanceError(f"{name} hash mismatch")
        if expected_hashes is not None and expected_hashes.get(name) != actual:
            raise ProvenanceError(f"configured {name} hash mismatch")
    if expected_hashes is not None:
        required = set(paths)
        if set(expected_hashes) != required:
            raise ProvenanceError("configured expected hashes must name every byte artifact")

    try:
        with np.load(mapping_path, allow_pickle=False) as artifact:
            full_axis = np.asarray(artifact["full_frequency_mhz"], dtype=np.float64)
            compact_axis = np.asarray(artifact["compact_frequency_mhz"], dtype=np.float64)
            mapping = np.asarray(artifact["full_to_compact"], dtype=np.int64)
        source_valid = np.asarray(np.load(source_valid_path, allow_pickle=False), dtype=bool)
        effective_mask = np.asarray(np.load(effective_mask_path, allow_pickle=False), dtype=bool)
    except (OSError, KeyError, ValueError) as exc:
        raise ProvenanceError(f"cannot load mapping inputs: {exc}") from exc
    if recorded_hashes.get("full_frequency_axis") != _array_sha256(full_axis):
        raise ProvenanceError("full frequency-axis hash mismatch")
    if recorded_hashes.get("compact_frequency_axis") != _array_sha256(compact_axis):
        raise ProvenanceError("compact frequency-axis hash mismatch")
    _validate_mapping_arrays(
        full_axis, compact_axis, source_valid, effective_mask, mapping
    )
    if not np.array_equal(np.asarray(spectrum.frequencies, dtype=np.float64), compact_axis):
        raise ProvenanceError("loaded compact frequency axis/order mismatch")
    if spectrum.power.shape[0] != compact_axis.size:
        raise ProvenanceError("loaded compact input row count mismatch")

    full_rows = np.flatnonzero(mapping >= 0)
    compact_mask = np.empty(compact_axis.size, dtype=bool)
    compact_mask[mapping[full_rows]] = effective_mask[full_rows]
    existing = np.ma.getmaskarray(spectrum.power)
    spectrum.power = np.ma.MaskedArray(
        spectrum.power.data,
        mask=existing | compact_mask[:, np.newaxis],
        copy=False,
    )
    return spectrum


def apply_configured_effective_mask(
    spectrum: DynamicSpectrum, config: dict
) -> DynamicSpectrum:
    """Apply configured ACF mask; an enabled required gate fails closed."""
    cfg = config.get("analysis", {}).get("bad_channel_mask")
    if not cfg:
        return spectrum
    if not cfg.get("enable", True):
        if cfg.get("required", False):
            raise ProvenanceError("required bad-channel mask is disabled")
        return spectrum
    required_fields = (
        "source_valid_path",
        "map_path",
        "effective_mask_path",
        "mapping_path",
        "provenance_path",
        "event",
        "instrument",
        "expected_hashes",
    )
    missing = [field for field in required_fields if not cfg.get(field)]
    if missing:
        if cfg.get("required", False):
            raise ProvenanceError(
                "required bad-channel provenance missing: " + ", ".join(missing)
            )
        return spectrum
    return apply_verified_effective_mask(
        spectrum,
        input_path=config["input_data_path"],
        source_valid_path=cfg["source_valid_path"],
        map_path=cfg["map_path"],
        effective_mask_path=cfg["effective_mask_path"],
        mapping_path=cfg["mapping_path"],
        provenance_path=cfg["provenance_path"],
        event=cfg["event"],
        instrument=cfg["instrument"],
        expected_hashes=cfg["expected_hashes"],
    )


def validate_configured_effective_mask(config: dict) -> None:
    """Revalidate configured artifacts before accepting a cached product."""
    cfg = config.get("analysis", {}).get("bad_channel_mask")
    if not cfg:
        return
    if not cfg.get("enable", True):
        if cfg.get("required", False):
            raise ProvenanceError("required bad-channel mask is disabled")
        return
    try:
        with np.load(config["input_data_path"], allow_pickle=False) as input_data:
            frequencies = np.asarray(input_data["frequencies_mhz"], dtype=float)
    except (OSError, KeyError, ValueError) as exc:
        raise ProvenanceError(f"cannot validate configured input axis: {exc}") from exc
    dummy = DynamicSpectrum(
        np.zeros((frequencies.size, 1), dtype=np.float32),
        frequencies,
        np.zeros(1),
    )
    apply_configured_effective_mask(dummy, config)


def configured_mask_cache_identity(config: dict) -> dict:
    """Return current artifact-byte identities for cache fingerprints."""
    cfg = config.get("analysis", {}).get("bad_channel_mask")
    if not cfg:
        return {"configured": False}
    identity = {
        "configured": True,
        "required": bool(cfg.get("required", False)),
        "event": cfg.get("event"),
        "instrument": cfg.get("instrument"),
    }
    for field in (
        "source_valid_path",
        "map_path",
        "effective_mask_path",
        "mapping_path",
        "provenance_path",
    ):
        path = cfg.get(field)
        try:
            identity[field] = _file_sha256(path) if path else None
        except OSError:
            identity[field] = "missing"
    return identity
