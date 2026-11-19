"""Artifact checks for the original Sze--Lee Chapter 3 figure atlas."""

from semphy.junctions import FIGURE_DIR


ATLAS_FILES = (
    "junction-atlas-formation.pdf",
    "junction-atlas-doping-bias.pdf",
    "junction-atlas-graded.pdf",
    "junction-atlas-capacitance-construction.pdf",
    "junction-atlas-injection.pdf",
    "junction-atlas-diode-characteristics.pdf",
    "junction-atlas-temperature-characteristics.pdf",
    "junction-atlas-small-signal.pdf",
    "junction-atlas-switching.pdf",
    "junction-atlas-breakdown-mechanisms.pdf",
    "junction-atlas-breakdown-design.pdf",
    "junction-atlas-curvature.pdf",
    "junction-atlas-heterojunction.pdf",
)


def test_chapter_three_atlas_contains_every_figure_plate() -> None:
    """Every original schematic plate referenced by Chapters 16--20 exists."""
    missing = [name for name in ATLAS_FILES if not (FIGURE_DIR / name).is_file()]
    assert missing == []
    assert all((FIGURE_DIR / name).stat().st_size > 5_000 for name in ATLAS_FILES)
