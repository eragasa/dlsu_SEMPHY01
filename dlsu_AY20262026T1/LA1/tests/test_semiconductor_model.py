"""Validation-test skeleton for the equilibrium semiconductor model."""

import numpy as np
import pytest

from src.semiconductor_model import BandValley, fermi_dirac


def test_fermi_dirac_is_stable_and_bounded() -> None:
    occupation = fermi_dirac(
        np.array([-1.0e6, 0.0, 1.0e6]),
        chemical_potential_ev=0.0,
        thermal_energy_ev=0.025,
    )
    assert np.all(np.isfinite(occupation))
    assert np.all((0.0 <= occupation) & (occupation <= 1.0))
    assert occupation[0] == pytest.approx(1.0)
    assert occupation[1] == pytest.approx(0.5)
    assert occupation[2] == pytest.approx(0.0)


def test_band_valley_dos_mass() -> None:
    valley = BandValley(
        label="test",
        edge_energy_ev=1.0,
        principal_masses_m0=(1.0, 2.0, 4.0),
    )
    assert valley.dos_mass_m0 == pytest.approx(2.0)


@pytest.mark.skip(reason="Part II implementation skeleton")
def test_n_increases_monotonically_with_mu() -> None:
    pass


@pytest.mark.skip(reason="Part II implementation skeleton")
def test_p_decreases_monotonically_with_mu() -> None:
    pass


@pytest.mark.skip(reason="Part II implementation skeleton")
def test_intrinsic_neutrality_residual() -> None:
    pass


@pytest.mark.skip(reason="Part II implementation skeleton")
def test_rigid_energy_shift_invariance() -> None:
    pass


@pytest.mark.skip(reason="Part II implementation skeleton")
def test_isotropic_mass_tensor_matches_scalar_formula() -> None:
    pass


@pytest.mark.skip(reason="Part II implementation skeleton")
def test_anisotropic_tensor_matches_determinant_formula() -> None:
    pass


@pytest.mark.skip(reason="Part II implementation skeleton")
def test_effective_densities_scale_as_temperature_to_three_halves() -> None:
    pass


@pytest.mark.skip(reason="Part II implementation skeleton")
def test_energy_cutoff_and_quadrature_convergence() -> None:
    pass
