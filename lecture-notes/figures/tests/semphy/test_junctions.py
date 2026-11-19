"""Validation tests for the p--n-junction teaching models."""

import numpy as np
import pytest

from semphy.junctions import (
    EPSILON_0,
    abrupt_depletion_profile,
    anderson_band_offsets,
    avalanche_breakdown_voltage,
    avalanche_integral,
    built_in_potential,
    depletion_capacitance,
    extract_one_sided_doping,
    ideal_diode_current_density,
    minority_diffusion_fd,
    solve_poisson_boltzmann_junction,
    transient_minority_diffusion_fd,
)


def test_abrupt_depletion_preserves_charge_and_voltage() -> None:
    """Depletion charge is neutral and its field integrates to the built-in drop."""
    result = abrupt_depletion_profile(5e22, 1e22, 1e16, 300.0, 11.7 * EPSILON_0)
    assert 5e22 * result["p_width_m"] == pytest.approx(
        1e22 * result["n_width_m"], rel=1e-14
    )
    assert result["potential_v"][-1] == pytest.approx(result["built_in_v"], rel=3e-3)


def test_poisson_boltzmann_matches_bulk_potential_difference() -> None:
    """Nonlinear finite differences recover the equilibrium contact potential."""
    expected = built_in_potential(5e22, 1e22, 1e16, 300.0)
    result = solve_poisson_boltzmann_junction(
        5e22, 1e22, 1e16, 300.0, 11.7 * EPSILON_0, points=241
    )
    drop = result["potential_v"][-1] - result["potential_v"][0]
    assert drop == pytest.approx(expected, rel=1e-8)
    assert result["maximum_residual"] < 1e-8


def test_capacitance_decreases_and_cv_recovers_doping() -> None:
    """An abrupt one-sided C--V curve has the expected monotonicity and slope."""
    reverse = np.linspace(0.0, 10.0, 101)
    acceptor, donor = 1e24, 2e21
    permittivity = 11.7 * EPSILON_0
    capacitance = depletion_capacitance(
        reverse, acceptor, donor, 1e16, 300.0, permittivity
    )
    extracted = extract_one_sided_doping(reverse, capacitance, permittivity)
    expected = acceptor * donor / (acceptor + donor)
    assert np.all(np.diff(capacitance) < 0.0)
    assert extracted == pytest.approx(np.full_like(extracted, expected), rel=2e-11)


def test_minority_diffusion_matches_finite_base_solution() -> None:
    """The finite-difference profile and edge flux converge to analytic values."""
    length, diffusivity, lifetime = 120e-6, 1e-3, 2e-6
    result = minority_diffusion_fd(length, diffusivity, lifetime, 1.0, points=401)
    position = result["position_m"]
    diffusion_length = np.sqrt(diffusivity * lifetime)
    exact = np.sinh((length - position) / diffusion_length) / np.sinh(
        length / diffusion_length
    )
    exact_flux = diffusivity / diffusion_length / np.tanh(length / diffusion_length)
    assert np.linalg.norm(result["excess_m3"] - exact) / np.linalg.norm(exact) < 2e-6
    assert result["particle_flux_m2_s"] == pytest.approx(exact_flux, rel=0.006)


def test_ideal_diode_has_zero_equilibrium_current() -> None:
    """The exponential diode equation vanishes at zero applied voltage."""
    current = ideal_diode_current_density(np.array([-1.0, 0.0, 0.1]), 1e-9, 300.0)
    assert current[1] == 0.0
    assert current[0] < 0.0 < current[2]


def test_transient_storage_is_positive_and_monotone() -> None:
    """Implicit removal maintains positivity and monotonically decreases charge."""
    result = transient_minority_diffusion_fd(
        120e-6, 1e-3, 2e-6, 1e20, np.linspace(0.0, 8e-6, 81), points=121
    )
    assert np.min(result["excess_m3"]) >= 0.0
    assert np.all(np.diff(result["stored_charge_c_m2"]) <= 1e-20)


def test_avalanche_root_satisfies_integral_criterion() -> None:
    """Breakdown root makes the equal-coefficient ionization integral unity."""
    parameters = (1e21, 0.8, 11.7 * EPSILON_0, 7e7, 1.23e8)
    breakdown = avalanche_breakdown_voltage(*parameters)
    assert avalanche_integral(breakdown, *parameters) == pytest.approx(1.0, rel=1e-10)


def test_avalanche_voltage_falls_with_doping() -> None:
    """The local triangular-field model predicts lower breakdown at higher doping."""
    parameters = (0.8, 11.7 * EPSILON_0, 7e7, 1.23e8)
    low = avalanche_breakdown_voltage(5e20, *parameters)
    high = avalanche_breakdown_voltage(2e21, *parameters)
    assert high < low


def test_anderson_offsets_obey_band_gap_identity() -> None:
    """Conduction and valence offsets reproduce the band-gap discontinuity."""
    conduction, valence = anderson_band_offsets(4.07, 1.42, 3.50, 2.16)
    assert conduction - valence == pytest.approx(2.16 - 1.42, abs=1e-14)
