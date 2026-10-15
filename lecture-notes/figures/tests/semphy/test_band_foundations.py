from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pytest
from semphy import band_foundations


@pytest.mark.parametrize(
    ("recipe", "filename"),
    [
        (band_foundations.piab_spectrum, "foundation-piab-spectrum.pdf"),
        (
            band_foundations.boundary_conditions,
            "foundation-boundary-conditions.pdf",
        ),
        (band_foundations.bloch_wave, "foundation-bloch-wave.pdf"),
        (
            band_foundations.two_site_splitting,
            "foundation-two-site-splitting.pdf",
        ),
        (band_foundations.levels_to_band, "foundation-levels-to-band.pdf"),
        (band_foundations.dimerized_bands, "foundation-dimerized-bands.pdf"),
    ],
)
def test_recipe_writes_a_nonempty_pdf(
    recipe: Callable[[], None],
    filename: str,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(band_foundations, "FIGURE_DIR", tmp_path)

    recipe()

    output = tmp_path / filename
    assert output.is_file()
    assert output.stat().st_size > 1_000
