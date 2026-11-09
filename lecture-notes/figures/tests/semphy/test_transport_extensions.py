"""Tests for the reduced semiconductor-transport extension models."""

import numpy as np
import pytest

from semphy.transport_extensions import (
    EPSILON_0,
    M_0,
    Q,
    excess_decay,
    hall_coefficient,
    ionization_rate,
    matthiessen_mobility,
    mott_gurney_current,
    quasi_fermi_splitting,
    rectangular_barrier_transmission,
    saturated_space_charge_current,
    saturation_velocity_model,
)


def test_matthiessen_mobility_is_below_each_limit() -> None:
    total = matthiessen_mobility(np.array([1000.0, 1500.0]), 2000.0)
    assert np.all(total < np.array([1000.0, 1500.0]))
    assert np.all(total < 2000.0)


def test_hall_coefficient_has_carrier_sign() -> None:
    density = 2.0e22
    assert hall_coefficient(density, -1) == pytest.approx(-1.0 / (Q * density))
    assert hall_coefficient(density, 1) == pytest.approx(1.0 / (Q * density))


def test_excess_decay_defines_the_lifetime() -> None:
    value = excess_decay(3.0e-6, 4.0e20, 3.0e-6)
    assert value == pytest.approx(4.0e20 / np.e)


def test_quasi_fermi_splitting_vanishes_at_equilibrium() -> None:
    intrinsic = 1.0e16
    electron = 1.0e22
    hole = intrinsic**2 / electron
    splitting = quasi_fermi_splitting(electron, hole, intrinsic, 300.0)
    assert splitting == pytest.approx(0.0, abs=1.0e-14)


def test_rectangular_barrier_transmission_is_bounded_and_width_sensitive() -> None:
    energy = np.linspace(0.05, 0.75, 100) * Q
    narrow = rectangular_barrier_transmission(energy, 0.8 * Q, 0.6e-9, 0.26 * M_0)
    wide = rectangular_barrier_transmission(energy, 0.8 * Q, 1.2e-9, 0.26 * M_0)
    assert np.all((narrow > 0.0) & (narrow < 1.0))
    assert np.all(wide < narrow)


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
