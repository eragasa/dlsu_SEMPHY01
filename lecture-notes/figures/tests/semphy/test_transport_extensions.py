"""Tests for the reduced semiconductor-transport extension models."""

import numpy as np
import pytest

from semphy.transport_extensions import (
    EPSILON_0,
    M_0,
    Q,
    dynamic_center_kinetics,
    excess_decay,
    fd_rectangular_barrier_transmission,
    hall_coefficient,
    ionization_rate,
    matthiessen_mobility,
    mott_gurney_current,
    poisson_potential_fd,
    quasi_fermi_splitting,
    rectangular_barrier_transmission,
    saturated_space_charge_current,
    saturation_velocity_model,
    space_charge_field_fd,
    srh_recombination_rate,
    steady_srh_profile_fd,
    transient_continuity_fd,
    two_carrier_hall_coefficient,
)


def test_matthiessen_mobility_is_below_each_limit() -> None:
    total = matthiessen_mobility(np.array([1000.0, 1500.0]), 2000.0)
    assert np.all(total < np.array([1000.0, 1500.0]))
    assert np.all(total < 2000.0)


def test_hall_coefficient_has_carrier_sign() -> None:
    density = 2.0e22
    assert hall_coefficient(density, -1) == pytest.approx(-1.0 / (Q * density))
    assert hall_coefficient(density, 1) == pytest.approx(1.0 / (Q * density))


def test_two_carrier_hall_reduces_to_single_carrier_limits() -> None:
    electron_only = two_carrier_hall_coefficient(2.0e22, 0.0, 0.12, 0.05)
    hole_only = two_carrier_hall_coefficient(0.0, 3.0e22, 0.12, 0.05)
    assert electron_only == pytest.approx(-1.0 / (Q * 2.0e22))
    assert hole_only == pytest.approx(1.0 / (Q * 3.0e22))


def test_two_carrier_hall_sign_uses_mobility_squared_weighting() -> None:
    coefficient = two_carrier_hall_coefficient(
        electron_density_m3=4.0e22,
        hole_density_m3=2.0e22,
        electron_mobility=0.05,
        hole_mobility=0.20,
    )
    assert coefficient > 0.0
    cancellation = two_carrier_hall_coefficient(4.0e22, 1.0e22, 0.10, 0.20)
    assert cancellation == pytest.approx(0.0, abs=1.0e-20)


def test_excess_decay_defines_the_lifetime() -> None:
    value = excess_decay(3.0e-6, 4.0e20, 3.0e-6)
    assert value == pytest.approx(4.0e20 / np.e)


def test_quasi_fermi_splitting_vanishes_at_equilibrium() -> None:
    intrinsic = 1.0e16
    electron = 1.0e22
    hole = intrinsic**2 / electron
    splitting = quasi_fermi_splitting(electron, hole, intrinsic, 300.0)
    assert splitting == pytest.approx(0.0, abs=1.0e-14)


def test_srh_rate_obeys_detailed_balance_and_sign() -> None:
    intrinsic = 1.0e16
    equilibrium_electron = 1.0e21
    equilibrium_hole = intrinsic**2 / equilibrium_electron
    equilibrium = srh_recombination_rate(
        equilibrium_electron,
        equilibrium_hole,
        intrinsic,
        2.0e-6,
        3.0e-6,
        intrinsic,
        intrinsic,
    )
    excess = srh_recombination_rate(
        equilibrium_electron + 1.0e19,
        equilibrium_hole + 1.0e19,
        intrinsic,
        2.0e-6,
        3.0e-6,
        intrinsic,
        intrinsic,
    )
    assert equilibrium == pytest.approx(0.0, abs=1.0)
    assert excess > 0.0


def test_transient_continuity_conserves_particles_apart_from_lifetime() -> None:
    position = np.linspace(-5.0e-3, 5.0e-3, 201, endpoint=False)
    times = np.linspace(0.0, 8.0e-6, 17)
    initial = np.exp(-0.5 * (position / 0.5e-3) ** 2)
    lifetime = 20.0e-6
    solution = transient_continuity_fd(
        position,
        times,
        initial,
        diffusivity_m2_s=1.0e-3,
        drift_velocity_m_s=40.0,
        lifetime_s=lifetime,
    )
    numerical_ratio = np.sum(solution[-1]) / np.sum(solution[0])
    assert numerical_ratio == pytest.approx(np.exp(-times[-1] / lifetime), rel=3.0e-5)


def test_steady_srh_profile_satisfies_surface_boundary() -> None:
    position = np.linspace(0.0, 40.0e-6, 61)
    diffusivity = 1.0e-3
    surface_velocity = 50.0
    profile = steady_srh_profile_fd(
        position,
        generation_m3_s=0.0,
        diffusivity_m2_s=diffusivity,
        equilibrium_electron_m3=1.0e21,
        equilibrium_hole_m3=1.0e11,
        intrinsic_m3=1.0e16,
        electron_lifetime_s=1.0e8,
        hole_lifetime_s=1.0e8,
        electron_auxiliary_m3=1.0e16,
        hole_auxiliary_m3=1.0e16,
        surface_velocity_m_s=surface_velocity,
        right_excess_m3=1.0e18,
    )
    dx = position[1] - position[0]
    surface_gradient = (profile[1] - profile[0]) / dx
    assert diffusivity * surface_gradient == pytest.approx(
        surface_velocity * profile[0], rel=2.0e-6
    )
    assert profile[-1] == pytest.approx(1.0e18, rel=2.0e-6)


def test_dynamic_center_kinetics_preserves_charge_invariant() -> None:
    times = np.linspace(0.0, 3.0e-6, 121)
    center_density = 2.0e20
    result = dynamic_center_kinetics(
        times,
        equilibrium_electron_m3=1.0e21,
        equilibrium_hole_m3=1.0e11,
        center_density_m3=center_density,
        electron_capture_m3_s=5.0e-15,
        hole_capture_m3_s=2.0e-14,
        electron_auxiliary_m3=1.0e16,
        hole_auxiliary_m3=1.0e16,
        generation_m3_s=2.0e26,
        pulse_end_s=0.5e-6,
    )
    invariant = (
        result["hole_m3"]
        - result["electron_m3"]
        - center_density * result["occupation"]
    )
    assert np.ptp(invariant) / abs(invariant[0]) < 1.0e-9
    assert np.all((result["occupation"] >= 0.0) & (result["occupation"] <= 1.0))


def test_rectangular_barrier_transmission_is_bounded_and_width_sensitive() -> None:
    energy = np.linspace(0.05, 0.75, 100) * Q
    narrow = rectangular_barrier_transmission(energy, 0.8 * Q, 0.6e-9, 0.26 * M_0)
    wide = rectangular_barrier_transmission(energy, 0.8 * Q, 1.2e-9, 0.26 * M_0)
    assert np.all((narrow > 0.0) & (narrow < 1.0))
    assert np.all(wide < narrow)


def test_finite_difference_barrier_converges_to_exact_transmission() -> None:
    energy = np.linspace(0.1, 0.7, 20) * Q
    exact = rectangular_barrier_transmission(energy, 0.8 * Q, 1.0e-9, 0.26 * M_0)
    coarse = fd_rectangular_barrier_transmission(
        energy, 0.8 * Q, 1.0e-9, 0.26 * M_0, 30
    )
    fine = fd_rectangular_barrier_transmission(
        energy, 0.8 * Q, 1.0e-9, 0.26 * M_0, 120
    )
    coarse_error = np.linalg.norm(coarse - exact)
    fine_error = np.linalg.norm(fine - exact)
    assert fine_error < coarse_error / 10.0
    assert fine_error / np.linalg.norm(exact) < 1.0e-4


def test_poisson_finite_difference_matches_uniform_charge_solution() -> None:
    position = np.linspace(0.0, 2.0e-6, 101)
    charge_density = 3.0e3
    permittivity = 11.7 * EPSILON_0
    potential = poisson_potential_fd(position, charge_density, permittivity, 0.0, 0.0)
    exact = (
        charge_density
        * position
        * (position[-1] - position)
        / (2.0 * permittivity)
    )
    assert np.max(np.abs(potential - exact)) / np.max(exact) < 1.0e-10


def test_space_charge_finite_difference_recovers_analytic_currents() -> None:
    voltage = 2.0
    length = 2.0e-6
    permittivity = 11.7 * EPSILON_0
    _, field_mu, current_mu = space_charge_field_fd(
        voltage,
        length,
        permittivity,
        81,
        mobility_m2_v_s=0.12,
    )
    _, field_sat, current_sat = space_charge_field_fd(
        voltage,
        length,
        permittivity,
        81,
        saturation_velocity_m_s=1.0e5,
    )
    expected_mu = mott_gurney_current(voltage, permittivity, 0.12, length)
    expected_sat = saturated_space_charge_current(
        voltage, permittivity, 1.0e5, length
    )
    assert current_mu == pytest.approx(expected_mu, rel=1.0e-3)
    assert current_sat == pytest.approx(expected_sat, rel=1.0e-8)
    assert np.all(np.diff(field_mu) >= 0.0)
    assert np.all(np.diff(field_sat) >= 0.0)


def test_space_charge_current_scaling() -> None:
    voltage = np.array([1.0, 2.0])
    permittivity = 11.7 * EPSILON_0
    quadratic = mott_gurney_current(voltage, permittivity, 0.1, 1.0e-6)
    linear = saturated_space_charge_current(voltage, permittivity, 1.0e5, 1.0e-6)
    assert quadratic[1] / quadratic[0] == pytest.approx(4.0)
    assert linear[1] / linear[0] == pytest.approx(2.0)


def test_saturation_velocity_model_has_correct_limits() -> None:
    mobility = 0.1
    saturation = 1.0e5
    low_field = 1.0e2
    high_field = 1.0e10
    velocity = saturation_velocity_model(
        np.array([low_field, high_field]), mobility, saturation, 2.0
    )
    assert velocity[0] == pytest.approx(mobility * low_field, rel=1.0e-6)
    assert velocity[1] == pytest.approx(saturation, rel=1.0e-6)


def test_ionization_rate_increases_with_field() -> None:
    rates = ionization_rate(np.array([1.0e7, 2.0e7, 4.0e7]), 7.0e7, 1.2e8)
    assert np.all(np.diff(rates) > 0.0)
