"""Validation checks for equilibrium-to-device integration relations."""

import numpy as np
import pytest

from semphy.device_prediction import (
    built_in_from_carrier_states,
    built_in_from_chemical_potentials,
    depletion_width_from_ionized_dopants,
    generalized_einstein_factor,
    local_log_sensitivity,
)


def test_built_in_voltage_is_energy_reference_invariant() -> None:
    """A common rigid energy shift leaves chemical-potential voltage unchanged."""
    reference = built_in_from_chemical_potentials(0.42, -0.38)
    for shift in (-5.0, -1.0, 0.5, 3.0):
        assert built_in_from_chemical_potentials(0.42 + shift, -0.38 + shift) == pytest.approx(reference)


def test_chemical_and_carrier_built_in_relations_agree() -> None:
    """Nondegenerate chemical-potential and carrier-density forms are equivalent."""
    temperature = 300.0
    intrinsic = 1e16
    electron = 1e22
    hole = 1e23
    carrier_voltage = built_in_from_carrier_states(electron, hole, intrinsic, temperature)
    thermal_voltage = 8.617333262e-5 * temperature
    mu_n = thermal_voltage * np.log(electron / intrinsic)
    mu_p = -thermal_voltage * np.log(hole / intrinsic)
    assert built_in_from_chemical_potentials(mu_n, mu_p) == pytest.approx(carrier_voltage, rel=2e-10)


def test_generalized_einstein_recovers_boltzmann_limit() -> None:
    """The degeneracy correction approaches unity far below the band edge."""
    factor = generalized_einstein_factor(np.array([-10.0, -8.0, 0.0, 5.0]))
    assert factor[0] == pytest.approx(1.0, rel=2e-5)
    assert factor[1] == pytest.approx(1.0, rel=2e-4)
    assert np.all(np.diff(factor) > 0.0)


def test_depletion_width_uses_ionized_charge_and_preserves_neutrality() -> None:
    """The width partition satisfies ionized depletion-charge neutrality."""
    total, p_width, n_width = depletion_width_from_ionized_dopants(
        0.8, 0.0, 5e22, 1e22, 11.7 * 8.8541878128e-12
    )
    assert total == pytest.approx(p_width + n_width)
    assert 5e22 * p_width == pytest.approx(1e22 * n_width, rel=1e-14)


def test_log_sensitivity_recovers_power_law_exponent() -> None:
    """Finite-difference logarithmic sensitivity recovers a known square law."""
    change = 0.01
    sensitivity = local_log_sensitivity(4.0, 4.0 * (1.0 + change) ** 2, change)
    assert sensitivity == pytest.approx(2.0)
