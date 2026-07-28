"""Tests for tau consistency and attribution matrix."""

import math

import pandas as pd

from galaxies.foreground.attribution_matrix import (
    ARCHIVE_SNAPSHOT_PATH,
    build_attribution_matrix,
)
from galaxies.foreground.tau_consistency import (
    consistency_status,
    scale_tau_1ghz_ms,
)


def test_scale_tau_alpha4():
    tau = scale_tau_1ghz_ms(1.0, 600.0, alpha=4.0)
    assert tau > 1.0


def test_consistency_gate_matches_physics():
    # C = 2π τ Δν ∈ (0.628, 12.57) for consistent
    assert consistency_status(0.05, 0.01) == "consistent"
    assert consistency_status(0.1, 50.0) == "inconsistent"


@__import__("pytest").mark.skipif(not __import__("pathlib").Path("galaxies/foreground/data/frozen_census/bursts.csv").exists(), reason="Faber2026 manuscript fixtures moved to the analysis repository")
def test_attribution_matrix_twelve_rows():
    df = build_attribution_matrix()
    assert len(df) == 12
    assert set(df.nickname) == {
        "zach",
        "whitney",
        "oran",
        "isha",
        "wilhelm",
        "phineas",
        "freya",
        "hamilton",
        "mahi",
        "chromatica",
        "casey",
        "johndoeii",
    }


@__import__("pytest").mark.skipif(not __import__("pathlib").Path("galaxies/foreground/data/frozen_census/bursts.csv").exists(), reason="Faber2026 manuscript fixtures moved to the analysis repository")
def test_freya_inverse_dnu_flagged():
    df = build_attribution_matrix()
    freya = df[df.nickname == "freya"].iloc[0]
    assert freya.dnu_status == "inverse_scaling"
    assert "inverse_dnu_scaling" in freya.multi_screen_triggers


@__import__("pytest").mark.skipif(not __import__("pathlib").Path("galaxies/foreground/data/frozen_census/bursts.csv").exists(), reason="Faber2026 manuscript fixtures moved to the analysis repository")
def test_unmeasured_dnu_placeholder():
    df = build_attribution_matrix()
    isha = df[df.nickname == "isha"].iloc[0]
    assert str(isha.dnu_chime_mhz).startswith("N/A —")


@__import__("pytest").mark.skipif(not __import__("pathlib").Path("galaxies/foreground/data/frozen_census/bursts.csv").exists(), reason="Faber2026 manuscript fixtures moved to the analysis repository")
def test_joint_tau_loaded_for_casey():
    df = build_attribution_matrix()
    casey = df[df.nickname == "casey"].iloc[0]
    assert math.isfinite(casey.tau_joint_1ghz_ms)


@__import__("pytest").mark.skipif(not __import__("pathlib").Path("galaxies/foreground/data/frozen_census/bursts.csv").exists(), reason="Faber2026 manuscript fixtures moved to the analysis repository")
def test_dm_budget_wired_to_registry():
    df = build_attribution_matrix()
    wilhelm = df[df.nickname == "wilhelm"].iloc[0]
    assert str(wilhelm.dm_budget_verdict).startswith("DM:")
    assert "sightline_budget not wired" not in wilhelm.dm_budget_verdict


@__import__("pytest").mark.skipif(not __import__("pathlib").Path("galaxies/foreground/data/frozen_census/bursts.csv").exists(), reason="Faber2026 manuscript fixtures moved to the analysis repository")
def test_wilhelm_two_screen_coherence_computed():
    # Preserve what the historical artifact claimed without treating a fresh
    # build from drifted inputs as current science.
    df = pd.read_csv(ARCHIVE_SNAPSHOT_PATH)
    wilhelm = df[df.nickname == "wilhelm"].iloc[0]
    assert "kpc^2" in str(wilhelm.two_screen_coherence)


@__import__("pytest").mark.skipif(not __import__("pathlib").Path("galaxies/foreground/data/frozen_census/bursts.csv").exists(), reason="Faber2026 manuscript fixtures moved to the analysis repository")
def test_archived_attribution_matrix_snapshot_is_complete():
    """The preserved snapshot remains inspectable without claiming freshness."""
    archived = pd.read_csv(ARCHIVE_SNAPSHOT_PATH)
    assert len(archived) == 12
    assert archived.nickname.nunique() == 12
    assert {"tau_joint_1ghz_ms", "dnu_chime_mhz", "dnu_dsa_mhz"} <= set(archived)
