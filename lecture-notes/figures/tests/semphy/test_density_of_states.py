from __future__ import annotations

import numpy as np
from semphy.density_of_states import (
    ELECTRON_MASS_KG,
    density_of_states_mass,
    gaussian_density_of_states,
    histogram_density_of_states,
    integrated_states,
    parabolic_band_energy,
    parabolic_density_of_states,
    sampled_band_density_of_states,
    tight_binding_band_1d,
    uniform_bz_grid_1d,
)


def test_tight_binding_grid_and_band_cover_first_bz() -> None:
    lattice_constant = 0.5e-9
    wave_vector = uniform_bz_grid_1d(lattice_constant, 20_000)
    energy = tight_binding_band_1d(
        wave_vector,
        lattice_constant_m=lattice_constant,
        onsite_ev=0.3,
        hopping_ev=1.2,
    )

    assert wave_vector.min() > -np.pi / lattice_constant
    assert wave_vector.max() < np.pi / lattice_constant
    assert np.isclose(energy.min(), 0.3 - 2.4, atol=1.0e-6)
    assert np.isclose(energy.max(), 0.3 + 2.4, atol=1.0e-6)


def test_gaussian_sampled_band_dos_preserves_declared_state_count() -> None:
    lattice_constant = 0.5e-9
    wave_vector = uniform_bz_grid_1d(lattice_constant, 8_000)
    band = tight_binding_band_1d(
        wave_vector,
        lattice_constant_m=lattice_constant,
    )
    energy_grid = np.linspace(-2.5, 2.5, 4_001)
    dos = sampled_band_density_of_states(
        band,
        energy_grid,
        sigma_ev=0.025,
        band_degeneracies=2.0,
    )

    assert np.isclose(integrated_states(energy_grid, dos), 2.0, rtol=1.0e-10)


def test_histogram_dos_preserves_explicit_weights() -> None:
    energies = np.array([-0.5, -0.1, 0.2, 0.6])
    weights = np.array([0.2, 0.3, 0.7, 0.8])
    edges = np.linspace(-1.0, 1.0, 9)
    centers, dos = histogram_density_of_states(
        energies,
        edges,
        state_weights=weights,
    )

    assert centers.shape == dos.shape
    assert np.isclose(np.sum(dos * np.diff(edges)), np.sum(weights))


def test_parabolic_band_and_dos_have_expected_dimension_dependence() -> None:
    energy = np.array([0.01, 0.04, 0.09])
    dos_1d = parabolic_density_of_states(
        energy,
        edge_ev=0.0,
        effective_mass_kg=0.25 * ELECTRON_MASS_KG,
        dimension=1,
    )
    dos_2d = parabolic_density_of_states(
        energy,
        edge_ev=0.0,
        effective_mass_kg=0.25 * ELECTRON_MASS_KG,
        dimension=2,
    )
    dos_3d = parabolic_density_of_states(
        energy,
        edge_ev=0.0,
        effective_mass_kg=0.25 * ELECTRON_MASS_KG,
        dimension=3,
    )

    assert np.allclose(dos_1d / dos_1d[0], [1.0, 0.5, 1.0 / 3.0])
    assert np.allclose(dos_2d / dos_2d[0], [1.0, 1.0, 1.0])
    assert np.allclose(dos_3d / dos_3d[0], [1.0, 2.0, 3.0])


def test_anisotropic_dos_mass_preserves_ellipsoid_product() -> None:
    masses = ELECTRON_MASS_KG * np.array([0.19, 0.19, 0.916])
    mass_dos = density_of_states_mass(masses)

    assert np.isclose(mass_dos**3, np.prod(masses), rtol=1.0e-14)
    assert np.isclose(mass_dos / ELECTRON_MASS_KG, 0.3215, rtol=2.0e-3)


def test_parabolic_band_uses_each_principal_mass() -> None:
    vectors = np.array([[1.0e9, 0.0], [0.0, 1.0e9]])
    energy = parabolic_band_energy(
        vectors,
        edge_ev=1.0,
        principal_masses_kg=[0.2 * ELECTRON_MASS_KG, 0.8 * ELECTRON_MASS_KG],
    )

    assert energy[0] > energy[1] > 1.0


def test_direct_gaussian_weights_are_not_implicitly_normalized() -> None:
    energies = np.array([-0.2, 0.4])
    weights = np.array([1.5, 0.5])
    grid = np.linspace(-1.0, 1.2, 4_001)
    dos = gaussian_density_of_states(
        energies,
        grid,
        sigma_ev=0.03,
        state_weights=weights,
    )

    assert np.isclose(integrated_states(grid, dos), 2.0, rtol=1.0e-10)
