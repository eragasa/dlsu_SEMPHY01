"""Validation tests for the notebook-colocated Bloch2D implementation."""

import numpy as np
import pytest

from notebooks.bloch2d import (
    CheckerboardBloch2D,
    gaussian_dos_2d,
    stable_fermi,
)


def test_hamiltonian_is_hermitian_off_symmetry() -> None:
    model = CheckerboardBloch2D()
    matrix = model.hamiltonian(np.array([0.173, -0.219]) * model.zone_edge)
    assert np.allclose(matrix, matrix.T.conj(), atol=1.0e-14)


def test_plane_wave_basis_size_matches_cutoff() -> None:
    model = CheckerboardBloch2D(reciprocal_mode_limit=3)
    matrix = model.hamiltonian([0.0, 0.0])
    assert matrix.shape == (49, 49)


def test_uniform_zone_is_endpoint_excluded_and_unique() -> None:
    model = CheckerboardBloch2D()
    axis, vectors = model.uniform_zone(12)
    assert len(vectors) == 144
    assert axis[0] == pytest.approx(-model.zone_edge)
    assert axis[-1] < model.zone_edge
    assert len(np.unique(vectors, axis=0)) == len(vectors)


def test_gaussian_dos_preserves_state_count() -> None:
    energies = np.array([0.0, 0.5, 1.0])
    grid = np.linspace(-0.5, 1.5, 5000)
    state_weight = 0.2
    density = gaussian_dos_2d(energies, grid, 0.03, state_weight)
    assert np.trapezoid(density, grid) == pytest.approx(
        len(energies) * state_weight,
        rel=1.0e-8,
    )


def test_fermi_occupation_is_stable_and_monotone() -> None:
    energy = np.array([-1.0e6, -0.1, 0.0, 0.1, 1.0e6])
    occupation = stable_fermi(energy, 0.0, 300.0)
    assert np.all(np.isfinite(occupation))
    assert occupation[0] == pytest.approx(1.0)
    assert occupation[-1] == pytest.approx(0.0)
    assert np.all(np.diff(occupation) <= 0.0)
