from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pytest
from semphy import density_of_states_figures


@pytest.mark.parametrize(
    ("recipe", "filename"),
    [
        (density_of_states_figures.band_to_dos, "dos-band-to-density.pdf"),
        (
            density_of_states_figures.dimensionality_comparison,
            "dos-dimensionality.pdf",
        ),
        (
            density_of_states_figures.numerical_convergence,
            "dos-numerical-convergence.pdf",
        ),
    ],
)
def test_recipe_writes_a_nonempty_pdf(
    recipe: Callable[[], None],
    filename: str,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(density_of_states_figures, "FIGURE_DIR", tmp_path)

    recipe()

    output = tmp_path / filename
    assert output.is_file()
    assert output.stat().st_size > 1_000
