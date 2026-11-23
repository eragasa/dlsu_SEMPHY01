"""Transparent silicon and GaAs examples for SEMPHY01 Long Assessment 2.

The parameter sets in this module are instructional controls, not universal
material truth. Students must replace them with values appropriate to the exact
material state inherited from Long Assessment 1. All calculations use SI units
unless a name explicitly contains ``cm3`` or ``cm2``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.constants import Boltzmann as K_B
from scipy.constants import elementary_charge as Q
from scipy.constants import epsilon_0 as EPSILON_0
from scipy.sparse import diags, eye
from scipy.sparse.linalg import spsolve


@dataclass(frozen=True)
class MaterialExample:
    """Store one readable instructional material parameter set.

    Parameters
    ----------
    name : str
        Material label.
    relative_permittivity : float
        Static relative dielectric constant.
    gap_zero_eV : float
        Varshni zero-temperature band gap in electron-volts.
    varshni_alpha_eV_K : float
        Varshni coefficient in electron-volts per kelvin.
    varshni_beta_K : float
        Varshni temperature in kelvin.
    conduction_dos_300_m3 : float
        Effective conduction density of states at 300 K in ``m^-3``.
    valence_dos_300_m3 : float
        Effective valence density of states at 300 K in ``m^-3``.
    electron_mobility_parameters : tuple[float, float, float, float]
        Caughey--Thomas ``(mu_min, mu_max, N_ref, alpha)`` with mobilities in
        ``cm^2 V^-1 s^-1`` and concentration in ``cm^-3``.
    hole_mobility_parameters : tuple[float, float, float, float]
        Hole Caughey--Thomas parameters in the same units.
    electron_lifetime_s : float
        Illustrative low-injection minority-electron lifetime in seconds.
    hole_lifetime_s : float
        Illustrative low-injection minority-hole lifetime in seconds.
    """

    name: str
    relative_permittivity: float
    gap_zero_eV: float
    varshni_alpha_eV_K: float
    varshni_beta_K: float
    conduction_dos_300_m3: float
    valence_dos_300_m3: float
    electron_mobility_parameters: tuple[float, float, float, float]
    hole_mobility_parameters: tuple[float, float, float, float]
    electron_lifetime_s: float
    hole_lifetime_s: float


SILICON = MaterialExample(
    name="Si",
    relative_permittivity=11.7,
    gap_zero_eV=1.170,
    varshni_alpha_eV_K=4.73e-4,
    varshni_beta_K=636.0,
    conduction_dos_300_m3=2.80e25,
    valence_dos_300_m3=1.04e25,
    electron_mobility_parameters=(52.2, 1417.0, 9.68e16, 0.680),
    hole_mobility_parameters=(44.9, 470.5, 2.23e17, 0.719),
    electron_lifetime_s=10.0e-6,
    hole_lifetime_s=10.0e-6,
)

GAAS = MaterialExample(
    name="GaAs",
    relative_permittivity=12.9,
    gap_zero_eV=1.519,
    varshni_alpha_eV_K=5.405e-4,
    varshni_beta_K=204.0,
    conduction_dos_300_m3=4.70e23,
    valence_dos_300_m3=7.00e24,
    electron_mobility_parameters=(800.0, 8500.0, 1.00e17, 0.70),
    hole_mobility_parameters=(50.0, 400.0, 2.00e17, 0.70),
    electron_lifetime_s=1.0e-9,
    hole_lifetime_s=1.0e-9,
)


def band_gap_eV(material: MaterialExample, temperature_K: np.ndarray | float) -> np.ndarray:
    """Evaluate the Varshni band-gap model.

    Parameters
    ----------
    material : MaterialExample
        Material parameter set.
    temperature_K : array_like
        Positive absolute temperature in kelvin.

    Returns
    -------
    numpy.ndarray
        Band gap in electron-volts.
    """
    temperature = np.asarray(temperature_K, dtype=float)
    if np.any(temperature <= 0.0):
        raise ValueError("temperature must be positive")
    return material.gap_zero_eV - (
        material.varshni_alpha_eV_K * temperature**2
        / (temperature + material.varshni_beta_K)
    )


def effective_densities_m3(
    material: MaterialExample, temperature_K: np.ndarray | float
) -> tuple[np.ndarray, np.ndarray]:
    """Scale effective densities of states with ``T^(3/2)``.

    Parameters
    ----------
    material : MaterialExample
        Material parameter set.
    temperature_K : array_like
        Positive absolute temperature in kelvin.

    Returns
    -------
    conduction_m3 : numpy.ndarray
        Effective conduction density in ``m^-3``.
    valence_m3 : numpy.ndarray
        Effective valence density in ``m^-3``.
    """
    temperature = np.asarray(temperature_K, dtype=float)
    if np.any(temperature <= 0.0):
        raise ValueError("temperature must be positive")
    scale = (temperature / 300.0) ** 1.5
    return material.conduction_dos_300_m3 * scale, material.valence_dos_300_m3 * scale


def intrinsic_concentration_m3(
    material: MaterialExample, temperature_K: np.ndarray | float
) -> np.ndarray:
    """Calculate nondegenerate intrinsic carrier concentration.

    Parameters
    ----------
    material : MaterialExample
        Material parameter set.
    temperature_K : array_like
        Positive absolute temperature in kelvin.

    Returns
    -------
    numpy.ndarray
        Intrinsic carrier concentration in ``m^-3``.
    """
    temperature = np.asarray(temperature_K, dtype=float)
    conduction, valence = effective_densities_m3(material, temperature)
    thermal_eV = K_B * temperature / Q
    return np.sqrt(conduction * valence) * np.exp(
        -band_gap_eV(material, temperature) / (2.0 * thermal_eV)
    )


def fully_ionized_equilibrium(
    material: MaterialExample,
    temperature_K: float,
    donor_m3: np.ndarray | float,
    acceptor_m3: np.ndarray | float,
) -> tuple[np.ndarray, np.ndarray]:
    """Solve nondegenerate mass action and fully ionized charge neutrality.

    Parameters
    ----------
    material : MaterialExample
        Material parameter set.
    temperature_K : float
        Positive absolute temperature in kelvin.
    donor_m3, acceptor_m3 : array_like
        Nonnegative donor and acceptor densities in ``m^-3``.

    Returns
    -------
    electron_m3 : numpy.ndarray
        Equilibrium electron density in ``m^-3``.
    hole_m3 : numpy.ndarray
        Equilibrium hole density in ``m^-3``.
    """
    donor = np.asarray(donor_m3, dtype=float)
    acceptor = np.asarray(acceptor_m3, dtype=float)
    if np.any(donor < 0.0) or np.any(acceptor < 0.0):
        raise ValueError("dopant densities must be nonnegative")
    net = donor - acceptor
    intrinsic = intrinsic_concentration_m3(material, temperature_K)
    radical = np.hypot(net, 2.0 * intrinsic)
    # Use algebraically equivalent roots that avoid subtracting nearly equal
    # numbers in strongly p-type or strongly n-type states.
    with np.errstate(divide="ignore", invalid="ignore"):
        electron = np.where(
            net >= 0.0,
            0.5 * (net + radical),
            2.0 * intrinsic**2 / (radical - net),
        )
        hole = np.where(
            net >= 0.0,
            2.0 * intrinsic**2 / (radical + net),
            0.5 * (-net + radical),
        )
    return electron, hole


def mobility_m2_Vs(
    material: MaterialExample,
    concentration_m3: np.ndarray | float,
    carrier: str,
) -> np.ndarray:
    """Evaluate a 300 K Caughey--Thomas mobility example.

    Parameters
    ----------
    material : MaterialExample
        Material parameter set.
    concentration_m3 : array_like
        Nonnegative total ionized-impurity concentration in ``m^-3``.
    carrier : {"electron", "hole"}
        Carrier whose drift mobility is requested.

    Returns
    -------
    numpy.ndarray
        Drift mobility in ``m^2 V^-1 s^-1``.
    """
    concentration = np.asarray(concentration_m3, dtype=float)
    if np.any(concentration < 0.0):
        raise ValueError("concentration must be nonnegative")
    if carrier == "electron":
        mu_min, mu_max, reference_cm3, exponent = material.electron_mobility_parameters
    elif carrier == "hole":
        mu_min, mu_max, reference_cm3, exponent = material.hole_mobility_parameters
    else:
        raise ValueError("carrier must be 'electron' or 'hole'")
    concentration_cm3 = concentration / 1.0e6
    mobility_cm2 = mu_min + (mu_max - mu_min) / (
        1.0 + (concentration_cm3 / reference_cm3) ** exponent
    )
    return mobility_cm2 * 1.0e-4


def conductivity_S_m(
    electron_m3: np.ndarray | float,
    hole_m3: np.ndarray | float,
    electron_mobility_m2_Vs: np.ndarray | float,
    hole_mobility_m2_Vs: np.ndarray | float,
) -> np.ndarray:
    """Calculate two-carrier conductivity.

    Parameters
    ----------
    electron_m3, hole_m3 : array_like
        Nonnegative carrier densities in ``m^-3``.
    electron_mobility_m2_Vs, hole_mobility_m2_Vs : array_like
        Nonnegative drift mobilities in ``m^2 V^-1 s^-1``.

    Returns
    -------
    numpy.ndarray
        Conductivity in ``S m^-1``.
    """
    return Q * (
        np.asarray(electron_m3) * np.asarray(electron_mobility_m2_Vs)
        + np.asarray(hole_m3) * np.asarray(hole_mobility_m2_Vs)
    )


def hall_coefficient_m3_C(
    electron_m3: np.ndarray | float,
    hole_m3: np.ndarray | float,
    electron_mobility_m2_Vs: np.ndarray | float,
    hole_mobility_m2_Vs: np.ndarray | float,
) -> np.ndarray:
    """Calculate the unit-scattering-factor two-carrier Hall coefficient.

    Parameters
    ----------
    electron_m3, hole_m3 : array_like
        Nonnegative carrier densities in ``m^-3``.
    electron_mobility_m2_Vs, hole_mobility_m2_Vs : array_like
        Nonnegative mobilities in ``m^2 V^-1 s^-1``.

    Returns
    -------
    numpy.ndarray
        Hall coefficient in ``m^3 C^-1``.
    """
    electron = np.asarray(electron_m3)
    hole = np.asarray(hole_m3)
    electron_mobility = np.asarray(electron_mobility_m2_Vs)
    hole_mobility = np.asarray(hole_mobility_m2_Vs)
    denominator = Q * (hole * hole_mobility + electron * electron_mobility) ** 2
    return (hole * hole_mobility**2 - electron * electron_mobility**2) / denominator


def diffusivity_m2_s(mobility_m2_Vs: np.ndarray | float, temperature_K: float) -> np.ndarray:
    """Apply the nondegenerate Einstein relation.

    Parameters
    ----------
    mobility_m2_Vs : array_like
        Mobility in ``m^2 V^-1 s^-1``.
    temperature_K : float
        Positive absolute temperature in kelvin.

    Returns
    -------
    numpy.ndarray
        Diffusivity in ``m^2 s^-1``.
    """
    if temperature_K <= 0.0:
        raise ValueError("temperature must be positive")
    return np.asarray(mobility_m2_Vs) * K_B * temperature_K / Q


def minority_profile_fd(
    length_m: float,
    diffusivity_m2_s_value: float,
    lifetime_s: float,
    boundary_excess_m3: float,
    *,
    points: int = 201,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Solve a steady finite-base minority diffusion problem.

    Parameters
    ----------
    length_m : float
        Positive neutral-region length in meters.
    diffusivity_m2_s_value : float
        Positive diffusivity in ``m^2 s^-1``.
    lifetime_s : float
        Positive lifetime in seconds.
    boundary_excess_m3 : float
        Nonnegative junction-edge excess density in ``m^-3``.
    points : int, optional
        Number of uniform grid nodes.

    Returns
    -------
    position_m : numpy.ndarray
        Spatial grid in meters.
    excess_m3 : numpy.ndarray
        Excess minority concentration in ``m^-3``.
    particle_flux_m2_s : float
        Junction-edge particle flux magnitude in ``m^-2 s^-1``.
    """
    if min(length_m, diffusivity_m2_s_value, lifetime_s) <= 0.0:
        raise ValueError("length, diffusivity, and lifetime must be positive")
    position = np.linspace(0.0, length_m, points)
    spacing = position[1] - position[0]
    interior = points - 2
    diffusion_length = np.sqrt(diffusivity_m2_s_value * lifetime_s)
    main = np.full(interior, -2.0 - spacing**2 / diffusion_length**2)
    matrix = diags((np.ones(interior - 1), main, np.ones(interior - 1)), (-1, 0, 1), format="csr")
    right_hand_side = np.zeros(interior)
    right_hand_side[0] = -boundary_excess_m3
    excess = np.zeros(points)
    excess[0] = boundary_excess_m3
    excess[1:-1] = spsolve(matrix, right_hand_side)
    particle_flux = diffusivity_m2_s_value * (excess[0] - excess[1]) / spacing
    return position, excess, float(particle_flux)


def transient_charge_removal_fd(
    length_m: float,
    diffusivity_m2_s_value: float,
    lifetime_s: float,
    initial_boundary_excess_m3: float,
    time_s: np.ndarray,
    *,
    points: int = 161,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Remove a steady minority profile with implicit Euler integration.

    Parameters
    ----------
    length_m, diffusivity_m2_s_value, lifetime_s : float
        Positive neutral-region length, diffusivity, and lifetime in SI units.
    initial_boundary_excess_m3 : float
        Nonnegative pre-switch boundary excess in ``m^-3``.
    time_s : numpy.ndarray
        Strictly increasing times beginning at zero.
    points : int, optional
        Number of spatial nodes.

    Returns
    -------
    position_m : numpy.ndarray
        Spatial grid in meters.
    excess_m3 : numpy.ndarray
        Excess-density history with shape ``(time, position)``.
    stored_charge_C_m2 : numpy.ndarray
        Stored charge per area in ``C m^-2``.
    """
    time = np.asarray(time_s, dtype=float)
    if time.ndim != 1 or time.size < 2 or time[0] != 0.0 or np.any(np.diff(time) <= 0.0):
        raise ValueError("time must increase strictly from zero")
    position, steady, _ = minority_profile_fd(
        length_m, diffusivity_m2_s_value, lifetime_s, initial_boundary_excess_m3, points=points
    )
    spacing = position[1] - position[0]
    interior = points - 2
    laplacian = diags((np.ones(interior - 1), -2.0 * np.ones(interior), np.ones(interior - 1)), (-1, 0, 1), format="csr") / spacing**2
    operator = diffusivity_m2_s_value * laplacian - eye(interior, format="csr") / lifetime_s
    history = np.zeros((time.size, points))
    history[0, 1:-1] = steady[1:-1]
    for index, step in enumerate(np.diff(time), start=1):
        history[index, 1:-1] = spsolve(
            eye(interior, format="csr") - step * operator,
            history[index - 1, 1:-1],
        )
    charge = Q * np.trapezoid(history, position, axis=1)
    return position, history, charge


def abrupt_junction(
    material: MaterialExample,
    temperature_K: float,
    acceptor_m3: float,
    donor_m3: float,
    applied_voltage_V: float = 0.0,
    *,
    points: int = 401,
) -> dict[str, np.ndarray | float]:
    """Calculate an abrupt fully ionized depletion profile.

    Parameters
    ----------
    material : MaterialExample
        Material parameter set.
    temperature_K : float
        Positive absolute temperature in kelvin.
    acceptor_m3, donor_m3 : float
        Positive ionized doping densities in ``m^-3``.
    applied_voltage_V : float, optional
        Junction voltage in volts, positive under forward bias.
    points : int, optional
        Number of depletion-region samples.

    Returns
    -------
    dict
        Built-in voltage, widths, position, charge, field, potential, and
        capacitance per area in SI units.
    """
    intrinsic = float(intrinsic_concentration_m3(material, temperature_K))
    thermal_voltage = K_B * temperature_K / Q
    built_in = thermal_voltage * np.log(acceptor_m3 * donor_m3 / intrinsic**2)
    drop = built_in - applied_voltage_V
    if drop <= 0.0:
        raise ValueError("forward voltage must remain below built-in voltage")
    permittivity = material.relative_permittivity * EPSILON_0
    width = np.sqrt(2.0 * permittivity * drop / Q * (1.0 / acceptor_m3 + 1.0 / donor_m3))
    p_width = donor_m3 / (acceptor_m3 + donor_m3) * width
    n_width = acceptor_m3 / (acceptor_m3 + donor_m3) * width
    position = np.linspace(-p_width, n_width, points)
    charge = np.where(position < 0.0, -Q * acceptor_m3, Q * donor_m3)
    field = np.where(
        position < 0.0,
        -Q * acceptor_m3 * (position + p_width) / permittivity,
        -Q * donor_m3 * (n_width - position) / permittivity,
    )
    potential = np.zeros_like(position)
    potential[1:] = np.cumsum(-0.5 * (field[1:] + field[:-1]) * np.diff(position))
    return {
        "built_in_V": built_in,
        "width_m": width,
        "p_width_m": p_width,
        "n_width_m": n_width,
        "position_m": position,
        "charge_C_m3": charge,
        "field_V_m": field,
        "potential_V": potential,
        "capacitance_F_m2": permittivity / width,
    }


def finite_base_saturation_current_A_m2(
    intrinsic_m3: float,
    acceptor_m3: float,
    donor_m3: float,
    electron_diffusivity_m2_s: float,
    hole_diffusivity_m2_s: float,
    electron_lifetime_s: float,
    hole_lifetime_s: float,
    p_neutral_length_m: float,
    n_neutral_length_m: float,
) -> float:
    """Calculate finite-base ideal-diode saturation-current density.

    Parameters
    ----------
    intrinsic_m3 : float
        Intrinsic carrier concentration in ``m^-3``.
    acceptor_m3, donor_m3 : float
        Positive ionized doping densities in ``m^-3``.
    electron_diffusivity_m2_s, hole_diffusivity_m2_s : float
        Positive minority diffusivities in ``m^2 s^-1``.
    electron_lifetime_s, hole_lifetime_s : float
        Positive minority lifetimes in seconds.
    p_neutral_length_m, n_neutral_length_m : float
        Positive p- and n-side neutral lengths in meters.

    Returns
    -------
    float
        Saturation-current density in ``A m^-2``.
    """
    electron_length = np.sqrt(electron_diffusivity_m2_s * electron_lifetime_s)
    hole_length = np.sqrt(hole_diffusivity_m2_s * hole_lifetime_s)
    electron_minority = intrinsic_m3**2 / acceptor_m3
    hole_minority = intrinsic_m3**2 / donor_m3
    return float(
        Q
        * (
            electron_diffusivity_m2_s
            * electron_minority
            / electron_length
            / np.tanh(p_neutral_length_m / electron_length)
            + hole_diffusivity_m2_s
            * hole_minority
            / hole_length
            / np.tanh(n_neutral_length_m / hole_length)
        )
    )


def ideal_diode_current_A_m2(
    voltage_V: np.ndarray | float,
    saturation_current_A_m2: float,
    temperature_K: float,
) -> np.ndarray:
    """Evaluate the ideal diode current density.

    Parameters
    ----------
    voltage_V : array_like
        Junction voltage in volts.
    saturation_current_A_m2 : float
        Positive saturation-current density in ``A m^-2``.
    temperature_K : float
        Positive absolute temperature in kelvin.

    Returns
    -------
    numpy.ndarray
        Current density in ``A m^-2``.
    """
    argument = np.clip(Q * np.asarray(voltage_V) / (K_B * temperature_K), -700.0, 700.0)
    return saturation_current_A_m2 * np.expm1(argument)
