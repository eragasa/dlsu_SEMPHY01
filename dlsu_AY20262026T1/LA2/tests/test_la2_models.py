"""Physical validation checks for the Si and GaAs LA2 examples."""

import numpy as np
import pytest

from la2_models import (
    GAAS,
    SILICON,
    abrupt_junction,
    conductivity_S_m,
    diffusivity_m2_s,
    finite_base_saturation_current_A_m2,
    fully_ionized_equilibrium,
    hall_coefficient_m3_C,
    intrinsic_concentration_m3,
    minority_profile_fd,
    mobility_m2_Vs,
    transient_charge_removal_fd,
)


def test_intrinsic_controls_have_expected_order() -> None:
    """The example intrinsic concentrations reproduce the Si/GaAs scale separation."""
    silicon = float(intrinsic_concentration_m3(SILICON, 300.0)) / 1e6
    gaas = float(intrinsic_concentration_m3(GAAS, 300.0)) / 1e6
    assert silicon == pytest.approx(6.12e9, rel=0.02)
    assert gaas == pytest.approx(2.04e6, rel=0.03)
    assert silicon > gaas


def test_mass_action_and_neutrality() -> None:
    """Fully ionized equilibrium satisfies mass action and charge neutrality."""
    concentration = np.logspace(20, 24, 20)
    intrinsic = float(intrinsic_concentration_m3(SILICON, 300.0))
    for donor, acceptor in ((concentration, 0.0), (0.0, concentration)):
        electron, hole = fully_ionized_equilibrium(
            SILICON, 300.0, donor, acceptor
        )
        np.testing.assert_allclose(electron * hole, intrinsic**2, rtol=1e-13)
        np.testing.assert_allclose(
            electron - hole, np.asarray(donor) - np.asarray(acceptor), rtol=1e-13
        )


def test_mobility_falls_with_ionized_impurity_density() -> None:
    """Both example Caughey--Thomas mobilities decrease with concentration."""
    concentration = np.logspace(20, 25, 100)
    for material in (SILICON, GAAS):
        for carrier in ("electron", "hole"):
            mobility = mobility_m2_Vs(material, concentration, carrier)
            assert np.all(np.diff(mobility) < 0.0)


def test_conductivity_and_hall_one_carrier_limit() -> None:
    """An n-type state has positive conductivity and negative Hall coefficient."""
    donor = 1e22
    electron, hole = fully_ionized_equilibrium(SILICON, 300.0, donor, 0.0)
    electron_mobility = mobility_m2_Vs(SILICON, donor, "electron")
    hole_mobility = mobility_m2_Vs(SILICON, donor, "hole")
    conductivity = conductivity_S_m(electron, hole, electron_mobility, hole_mobility)
    hall = hall_coefficient_m3_C(electron, hole, electron_mobility, hole_mobility)
    assert conductivity > 0.0
    assert hall < 0.0


def test_einstein_relation_is_dimensionally_consistent() -> None:
    """Diffusivity divided by mobility equals the thermal voltage."""
    mobility = np.array([0.01, 0.1, 1.0])
    diffusivity = diffusivity_m2_s(mobility, 300.0)
    np.testing.assert_allclose(diffusivity / mobility, 0.0258519998, rtol=1e-7)


def test_minority_profile_converges_to_finite_base_solution() -> None:
    """The finite-difference minority profile matches the analytic solution."""
    diffusivity, lifetime = 1e-3, 2e-6
    diffusion_length = np.sqrt(diffusivity * lifetime)
    length = 4 * diffusion_length
    position, excess, flux = minority_profile_fd(length, diffusivity, lifetime, 1.0, points=401)
    exact = np.sinh((length - position) / diffusion_length) / np.sinh(length / diffusion_length)
    exact_flux = diffusivity / diffusion_length / np.tanh(length / diffusion_length)
    assert np.linalg.norm(excess - exact) / np.linalg.norm(exact) < 3e-6
    assert flux == pytest.approx(exact_flux, rel=0.006)


def test_transient_removal_is_positive_and_monotone() -> None:
    """Implicit removal preserves positivity and decreases stored charge."""
    time = np.linspace(0.0, 12e-6, 121)
    _, history, charge = transient_charge_removal_fd(120e-6, 1e-3, 2e-6, 1e20, time)
    assert np.min(history) >= 0.0
    assert np.all(np.diff(charge) <= 1e-20)


def test_depletion_solution_preserves_charge_and_capacitance_identity() -> None:
    """Abrupt depletion widths are neutral and satisfy C times W equals epsilon."""
    result = abrupt_junction(SILICON, 300.0, 1e23, 1e22)
    assert 1e23 * result["p_width_m"] == pytest.approx(1e22 * result["n_width_m"], rel=1e-14)
    permittivity = SILICON.relative_permittivity * 8.8541878128e-12
    assert result["capacitance_F_m2"] * result["width_m"] == pytest.approx(permittivity)


def test_gaas_example_has_smaller_ideal_saturation_current() -> None:
    """The declared controls give lower ideal saturation current for GaAs."""
    currents = []
    for material in (SILICON, GAAS):
        intrinsic = float(intrinsic_concentration_m3(material, 300.0))
        electron_diffusivity = float(diffusivity_m2_s(mobility_m2_Vs(material, 1e23, "electron"), 300.0))
        hole_diffusivity = float(diffusivity_m2_s(mobility_m2_Vs(material, 1e22, "hole"), 300.0))
        electron_length = np.sqrt(electron_diffusivity * material.electron_lifetime_s)
        hole_length = np.sqrt(hole_diffusivity * material.hole_lifetime_s)
        currents.append(
            finite_base_saturation_current_A_m2(
                intrinsic,
                1e23,
                1e22,
                electron_diffusivity,
                hole_diffusivity,
                material.electron_lifetime_s,
                material.hole_lifetime_s,
                4 * electron_length,
                4 * hole_length,
            )
        )
    assert currents[1] < currents[0]
