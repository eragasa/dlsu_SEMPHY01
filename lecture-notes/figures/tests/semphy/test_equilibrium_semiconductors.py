from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pytest
from semphy import equilibrium_semiconductors


@pytest.mark.parametrize(
    ("recipe", "filename"),
    [
        (
            equilibrium_semiconductors.band_formation_schematic,
            "intrinsic-band-formation.pdf",
        ),
        (
            equilibrium_semiconductors.kspace_state_counting,
            "intrinsic-kspace-state-counting.pdf",
        ),
        (
            equilibrium_semiconductors.states_occupation_population,
            "intrinsic-states-occupation.pdf",
        ),
        (
            equilibrium_semiconductors.fermi_boltzmann_comparison,
            "intrinsic-fermi-boltzmann.pdf",
        ),
        (
            equilibrium_semiconductors.intrinsic_temperature_dependence,
            "intrinsic-temperature-dependence.pdf",
        ),
        (
            equilibrium_semiconductors.extrinsic_levels_and_charge_states,
            "extrinsic-dopant-levels.pdf",
        ),
        (
            equilibrium_semiconductors.extrinsic_ionization_fractions,
            "extrinsic-ionization-fractions.pdf",
        ),
        (
            equilibrium_semiconductors.extrinsic_compensation,
            "extrinsic-compensation.pdf",
        ),
        (
            equilibrium_semiconductors.extrinsic_temperature_regimes,
            "extrinsic-temperature-regimes.pdf",
        ),
    ],
)
def test_recipe_writes_a_nonempty_pdf(
    recipe: Callable[[], None],
    filename: str,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(equilibrium_semiconductors, "FIGURE_DIR", tmp_path)

    recipe()

    output = tmp_path / filename
    assert output.is_file()
    assert output.stat().st_size > 1_000
