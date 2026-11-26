"""Numerical models and deterministic figures for the Sze--Lee p--n junction.

The module supports Chapters 16--20 of the SEMPHY01 notes.  It combines
closed-form depletion results with finite-difference electrostatics, minority-
carrier diffusion, transient charge storage, capacitance profiling, avalanche
integration, and a simple heterojunction band-offset construction.

All public functions use SI units unless an argument name explicitly ends in
``_ev``.  The models are transparent teaching calculations rather than a
production semiconductor-device simulator.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq, root
from scipy.sparse import diags, eye
from scipy.sparse.linalg import spsolve

Q = 1.602_176_634e-19
K_B = 1.380_649e-23
EPSILON_0 = 8.854_187_8128e-12
FIGURE_DIR = Path(__file__).resolve().parents[2]

mpl.rcParams.update(
    {
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "legend.fontsize": 8,
        "figure.dpi": 140,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.04,
        "pdf.fonttype": 42,
    }
)


def _save(figure: plt.Figure, name: str) -> None:
    """Write a timestamp-free PDF beneath the manuscript figure directory."""
    metadata = {
        "Creator": "SEMPHY01 deterministic junction figure generator",
        "Producer": "Matplotlib",
        "CreationDate": None,
        "ModDate": None,
    }
    figure.savefig(FIGURE_DIR / name, metadata=metadata)
    plt.close(figure)


def built_in_potential(
    acceptor_m3: float,
    donor_m3: float,
    intrinsic_m3: float,
    temperature_k: float,
) -> float:
    """Calculate the nondegenerate abrupt-junction built-in potential.

    Parameters
    ----------
    acceptor_m3, donor_m3 : float
        Positive ionized acceptor and donor densities in ``m^-3``.
    intrinsic_m3 : float
        Positive intrinsic carrier density in ``m^-3``.
    temperature_k : float
        Positive absolute temperature in kelvin.

    Returns
    -------
    float
        Built-in potential in volts.

    Raises
    ------
    ValueError
        If a density or temperature is nonpositive.
    """
    if min(acceptor_m3, donor_m3, intrinsic_m3, temperature_k) <= 0.0:
        raise ValueError("densities and temperature must be positive")
    thermal_voltage = K_B * temperature_k / Q
    return thermal_voltage * np.log(
        acceptor_m3 * donor_m3 / intrinsic_m3**2
    )


def abrupt_depletion_profile(
    acceptor_m3: float,
    donor_m3: float,
    intrinsic_m3: float,
    temperature_k: float,
    permittivity_f_m: float,
    applied_voltage_v: float = 0.0,
    *,
    points: int = 401,
) -> dict[str, np.ndarray | float]:
    """Return the depletion-approximation solution of an abrupt junction.

    Parameters
    ----------
    acceptor_m3, donor_m3 : float
        Positive abrupt doping densities in ``m^-3``.
    intrinsic_m3 : float
        Positive intrinsic density in ``m^-3``.
    temperature_k : float
        Positive absolute temperature in kelvin.
    permittivity_f_m : float
        Positive semiconductor permittivity in ``F m^-1``.
    applied_voltage_v : float, optional
        Junction voltage, positive for forward bias and negative for reverse
        bias, in volts.
    points : int, optional
        Number of spatial samples across the depletion region.

    Returns
    -------
    dict
        Position, charge density, electric field, electrostatic potential,
        p-side width, n-side width, total width, maximum-field magnitude, and
        built-in potential.

    Raises
    ------
    ValueError
        If parameters are invalid or forward bias removes the depletion drop.

    Notes
    -----
    Charge neutrality imposes ``N_A*x_p=N_D*x_n``.  The depletion voltage is
    ``V_bi-V`` under the voltage convention used here.
    """
    if min(acceptor_m3, donor_m3, intrinsic_m3, temperature_k, permittivity_f_m) <= 0.0:
        raise ValueError("densities, temperature, and permittivity must be positive")
    if type(points) is not int or points < 20:
        raise ValueError("points must be an integer of at least 20")
    built_in = built_in_potential(
        acceptor_m3, donor_m3, intrinsic_m3, temperature_k
    )
    depletion_voltage = built_in - applied_voltage_v
    if depletion_voltage <= 0.0:
        raise ValueError("applied forward voltage must remain below built-in potential")

    total_width = np.sqrt(
        2.0
        * permittivity_f_m
        * depletion_voltage
        / Q
        * (1.0 / acceptor_m3 + 1.0 / donor_m3)
    )
    p_width = donor_m3 / (acceptor_m3 + donor_m3) * total_width
    n_width = acceptor_m3 / (acceptor_m3 + donor_m3) * total_width
    position = np.linspace(-p_width, n_width, points)
    charge = np.where(position < 0.0, -Q * acceptor_m3, Q * donor_m3)
    field = np.where(
        position < 0.0,
        -Q * acceptor_m3 * (position + p_width) / permittivity_f_m,
        -Q * donor_m3 * (n_width - position) / permittivity_f_m,
    )
    # Integrate -E=dpsi/dx and normalize the p-side potential to zero.
    potential = cumulative_trapezoid(-field, position, initial=0.0)
    maximum_field = Q * acceptor_m3 * p_width / permittivity_f_m
    return {
        "position_m": position,
        "charge_c_m3": charge,
        "field_v_m": field,
        "potential_v": potential,
        "p_width_m": p_width,
        "n_width_m": n_width,
        "width_m": total_width,
        "maximum_field_v_m": maximum_field,
        "built_in_v": built_in,
    }


def solve_poisson_boltzmann_junction(
    acceptor_m3: float,
    donor_m3: float,
    intrinsic_m3: float,
    temperature_k: float,
    permittivity_f_m: float,
    *,
    points: int = 301,
    extent_factor: float = 4.0,
) -> dict[str, np.ndarray | float]:
    """Solve equilibrium junction electrostatics without the depletion approximation.

    Parameters
    ----------
    acceptor_m3, donor_m3 : float
        Positive abrupt doping densities in ``m^-3``.
    intrinsic_m3 : float
        Positive intrinsic density in ``m^-3``.
    temperature_k : float
        Positive absolute temperature in kelvin.
    permittivity_f_m : float
        Positive semiconductor permittivity in ``F m^-1``.
    points : int, optional
        Number of uniform finite-difference nodes.
    extent_factor : float, optional
        Domain length on each side in multiples of the corresponding depletion
        width.

    Returns
    -------
    dict
        Position, electrostatic potential, field, mobile densities, charge
        density, and nonlinear residual.

    Raises
    ------
    ValueError
        If parameters or grid controls are invalid.
    RuntimeError
        If the nonlinear finite-difference solve fails.

    Notes
    -----
    Boltzmann carriers obey ``n=n_i exp(psi/V_T)`` and
    ``p=n_i exp(-psi/V_T)``.  The neutral bulk potentials are obtained from
    ``n-p=N_D-N_A``.  Interior rows discretize Poisson's equation.
    """
    if type(points) is not int or points < 40:
        raise ValueError("points must be an integer of at least 40")
    if extent_factor <= 1.0:
        raise ValueError("extent_factor must exceed one")
    depletion = abrupt_depletion_profile(
        acceptor_m3,
        donor_m3,
        intrinsic_m3,
        temperature_k,
        permittivity_f_m,
    )
    thermal_voltage = K_B * temperature_k / Q
    left = -extent_factor * float(depletion["p_width_m"])
    right = extent_factor * float(depletion["n_width_m"])
    position = np.linspace(left, right, points)
    dx = position[1] - position[0]
    net_doping = np.where(position < 0.0, -acceptor_m3, donor_m3)
    left_potential = thermal_voltage * np.arcsinh(-acceptor_m3 / (2.0 * intrinsic_m3))
    right_potential = thermal_voltage * np.arcsinh(donor_m3 / (2.0 * intrinsic_m3))
    p_width = float(depletion["p_width_m"])
    n_width = float(depletion["n_width_m"])
    initial = np.empty(points)
    initial[position <= -p_width] = left_potential
    initial[position >= n_width] = right_potential
    p_mask = (position > -p_width) & (position < 0.0)
    n_mask = (position >= 0.0) & (position < n_width)
    initial[p_mask] = left_potential + (
        Q * acceptor_m3 * (position[p_mask] + p_width) ** 2
        / (2.0 * permittivity_f_m)
    )
    depletion_drop = float(depletion["built_in_v"])
    initial[n_mask] = left_potential + depletion_drop - (
        Q * donor_m3 * (n_width - position[n_mask]) ** 2
        / (2.0 * permittivity_f_m)
    )
    # Match the exact neutral-boundary voltage while preserving the depletion shape.
    initial = left_potential + (initial - left_potential) * (
        (right_potential - left_potential) / depletion_drop
    )

    def residual(potential: np.ndarray) -> np.ndarray:
        reduced = np.clip(potential / thermal_voltage, -80.0, 80.0)
        electron = intrinsic_m3 * np.exp(reduced)
        hole = intrinsic_m3 * np.exp(-reduced)
        charge = Q * (hole - electron + net_doping)
        value = np.empty_like(potential)
        value[0] = (potential[0] - left_potential) / thermal_voltage
        value[-1] = (potential[-1] - right_potential) / thermal_voltage
        value[1:-1] = (
            potential[:-2]
            - 2.0 * potential[1:-1]
            + potential[2:]
            + dx**2 * charge[1:-1] / permittivity_f_m
        ) / thermal_voltage
        return value

    solution = root(residual, initial, method="hybr", options={"xtol": 1.0e-9, "maxfev": 10000})
    maximum_residual = float(np.max(np.abs(residual(solution.x))))
    if not solution.success or maximum_residual > 2.0e-7:
        raise RuntimeError("Poisson--Boltzmann junction solve did not converge")
    potential = solution.x
    reduced = np.clip(potential / thermal_voltage, -80.0, 80.0)
    electron = intrinsic_m3 * np.exp(reduced)
    hole = intrinsic_m3 * np.exp(-reduced)
    charge = Q * (hole - electron + net_doping)
    field = -np.gradient(potential, position, edge_order=2)
    return {
        "position_m": position,
        "potential_v": potential,
        "field_v_m": field,
        "electron_m3": electron,
        "hole_m3": hole,
        "charge_c_m3": charge,
        "maximum_residual": maximum_residual,
    }


def depletion_capacitance(
    reverse_voltage_v: np.ndarray | float,
    acceptor_m3: float,
    donor_m3: float,
    intrinsic_m3: float,
    temperature_k: float,
    permittivity_f_m: float,
) -> np.ndarray:
    """Calculate abrupt-junction depletion capacitance per unit area.

    Parameters
    ----------
    reverse_voltage_v : array_like
        Nonnegative reverse-bias magnitude in volts.
    acceptor_m3, donor_m3, intrinsic_m3 : float
        Positive carrier densities in ``m^-3``.
    temperature_k : float
        Positive absolute temperature in kelvin.
    permittivity_f_m : float
        Positive permittivity in ``F m^-1``.

    Returns
    -------
    numpy.ndarray
        Depletion capacitance per area in ``F m^-2``.

    Raises
    ------
    ValueError
        If reverse bias is negative or a parameter is nonpositive.
    """
    reverse = np.asarray(reverse_voltage_v, dtype=float)
    if np.any(reverse < 0.0):
        raise ValueError("reverse voltage must be nonnegative")
    if permittivity_f_m <= 0.0:
        raise ValueError("permittivity must be positive")
    built_in = built_in_potential(acceptor_m3, donor_m3, intrinsic_m3, temperature_k)
    width = np.sqrt(
        2.0
        * permittivity_f_m
        * (built_in + reverse)
        / Q
        * (1.0 / acceptor_m3 + 1.0 / donor_m3)
    )
    return permittivity_f_m / width


def extract_one_sided_doping(
    reverse_voltage_v: np.ndarray,
    capacitance_f_m2: np.ndarray,
    permittivity_f_m: float,
) -> np.ndarray:
    """Extract one-sided junction doping from a numerical C--V curve.

    Parameters
    ----------
    reverse_voltage_v : numpy.ndarray
        Strictly increasing reverse-bias values in volts.
    capacitance_f_m2 : numpy.ndarray
        Positive capacitance per area in ``F m^-2``.
    permittivity_f_m : float
        Positive semiconductor permittivity in ``F m^-1``.

    Returns
    -------
    numpy.ndarray
        Extracted lightly doped concentration in ``m^-3``.

    Raises
    ------
    ValueError
        If arrays are incompatible or parameters are invalid.

    Notes
    -----
    The one-sided abrupt formula is
    ``N = 2/[q*epsilon*d(1/C**2)/dV_R]``.
    """
    voltage = np.asarray(reverse_voltage_v, dtype=float)
    capacitance = np.asarray(capacitance_f_m2, dtype=float)
    if voltage.ndim != 1 or voltage.size < 4 or capacitance.shape != voltage.shape:
        raise ValueError("voltage and capacitance must be matching one-dimensional arrays")
    if np.any(np.diff(voltage) <= 0.0) or np.any(capacitance <= 0.0):
        raise ValueError("voltage must increase and capacitance must be positive")
    if permittivity_f_m <= 0.0:
        raise ValueError("permittivity must be positive")
    slope = np.gradient(1.0 / capacitance**2, voltage, edge_order=2)
    return 2.0 / (Q * permittivity_f_m * slope)


def minority_diffusion_fd(
    length_m: float,
    diffusivity_m2_s: float,
    lifetime_s: float,
    junction_excess_m3: float,
    *,
    points: int = 201,
) -> dict[str, np.ndarray | float]:
    """Solve the steady neutral-region minority diffusion equation.

    Parameters
    ----------
    length_m : float
        Positive neutral-region length in meters.
    diffusivity_m2_s : float
        Positive minority diffusivity in ``m^2 s^-1``.
    lifetime_s : float
        Positive minority lifetime in seconds.
    junction_excess_m3 : float
        Nonnegative excess density at the depletion edge in ``m^-3``.
    points : int, optional
        Number of finite-difference nodes.

    Returns
    -------
    dict
        Position, excess-density profile, and junction particle flux in
        ``m^-2 s^-1``.

    Raises
    ------
    ValueError
        If parameters are invalid.
    """
    if min(length_m, diffusivity_m2_s, lifetime_s) <= 0.0 or junction_excess_m3 < 0.0:
        raise ValueError("length, diffusivity, and lifetime must be positive")
    if type(points) is not int or points < 4:
        raise ValueError("points must be an integer of at least four")
    position = np.linspace(0.0, length_m, points)
    dx = position[1] - position[0]
    interior = points - 2
    main = np.full(interior, -2.0 - dx**2 / (diffusivity_m2_s * lifetime_s))
    matrix = diags((np.ones(interior - 1), main, np.ones(interior - 1)), (-1, 0, 1), format="csr")
    right_hand_side = np.zeros(interior)
    right_hand_side[0] -= junction_excess_m3
    excess = np.empty(points)
    excess[0] = junction_excess_m3
    excess[-1] = 0.0
    excess[1:-1] = spsolve(matrix, right_hand_side)
    particle_flux = diffusivity_m2_s * (excess[0] - excess[1]) / dx
    return {"position_m": position, "excess_m3": excess, "particle_flux_m2_s": particle_flux}


def ideal_diffusion_saturation_current_density(
    acceptor_m3: float,
    donor_m3: float,
    intrinsic_m3: float,
    electron_diffusivity_m2_s: float,
    hole_diffusivity_m2_s: float,
    electron_lifetime_s: float,
    hole_lifetime_s: float,
) -> float:
    """Calculate the long-base diffusion saturation-current density.

    Parameters
    ----------
    acceptor_m3, donor_m3, intrinsic_m3 : float
        Positive acceptor, donor, and intrinsic carrier densities in ``m^-3``.
    electron_diffusivity_m2_s, hole_diffusivity_m2_s : float
        Positive minority-electron and minority-hole diffusivities in
        ``m^2 s^-1``.
    electron_lifetime_s, hole_lifetime_s : float
        Positive minority-carrier lifetimes in seconds.

    Returns
    -------
    float
        Saturation-current density in ``A m^-2``.

    Notes
    -----
    The long-base result is
    ``q*n_i^2*(D_n/(L_n*N_A) + D_p/(L_p*N_D))``, with
    ``L_n=sqrt(D_n*tau_n)`` and ``L_p=sqrt(D_p*tau_p)``.
    """
    parameters = (
        acceptor_m3,
        donor_m3,
        intrinsic_m3,
        electron_diffusivity_m2_s,
        hole_diffusivity_m2_s,
        electron_lifetime_s,
        hole_lifetime_s,
    )
    if min(parameters) <= 0.0:
        raise ValueError("densities, diffusivities, and lifetimes must be positive")
    electron_diffusion_length = np.sqrt(
        electron_diffusivity_m2_s * electron_lifetime_s
    )
    hole_diffusion_length = np.sqrt(hole_diffusivity_m2_s * hole_lifetime_s)
    return float(
        Q
        * intrinsic_m3**2
        * (
            electron_diffusivity_m2_s
            / (electron_diffusion_length * acceptor_m3)
            + hole_diffusivity_m2_s / (hole_diffusion_length * donor_m3)
        )
    )


def ideal_diode_current_density(
    voltage_v: np.ndarray | float,
    saturation_current_a_m2: float,
    temperature_k: float,
    *,
    ideality: float = 1.0,
) -> np.ndarray:
    """Calculate an exponential diode current density.

    Parameters
    ----------
    voltage_v : array_like
        Junction voltage in volts.
    saturation_current_a_m2 : float
        Positive saturation-current density in ``A m^-2``.
    temperature_k : float
        Positive absolute temperature in kelvin.
    ideality : float, optional
        Positive ideality factor.

    Returns
    -------
    numpy.ndarray
        Current density in ``A m^-2``.

    Raises
    ------
    ValueError
        If current, temperature, or ideality is nonpositive.
    """
    if saturation_current_a_m2 <= 0.0 or temperature_k <= 0.0 or ideality <= 0.0:
        raise ValueError("saturation current, temperature, and ideality must be positive")
    voltage = np.asarray(voltage_v, dtype=float)
    argument = np.clip(Q * voltage / (ideality * K_B * temperature_k), -700.0, 700.0)
    return saturation_current_a_m2 * np.expm1(argument)


def transient_minority_diffusion_fd(
    length_m: float,
    diffusivity_m2_s: float,
    lifetime_s: float,
    initial_junction_excess_m3: float,
    time_s: np.ndarray,
    *,
    points: int = 201,
) -> dict[str, np.ndarray]:
    """Integrate minority-carrier removal after a junction boundary is switched off.

    Parameters
    ----------
    length_m : float
        Positive neutral-region length in meters.
    diffusivity_m2_s : float
        Positive diffusivity in ``m^2 s^-1``.
    lifetime_s : float
        Positive lifetime in seconds.
    initial_junction_excess_m3 : float
        Nonnegative pre-switch excess density at the junction in ``m^-3``.
    time_s : numpy.ndarray
        Strictly increasing times in seconds beginning at zero.
    points : int, optional
        Number of spatial nodes.

    Returns
    -------
    dict
        Position, time-dependent excess profile, stored charge per area in
        ``C m^-2``, and extraction-current density in ``A m^-2``.

    Raises
    ------
    ValueError
        If physical or grid parameters are invalid.

    Notes
    -----
    The initial profile is the steady finite-base diffusion solution.  At
    ``t=0`` both boundaries are set to equilibrium.  Backward Euler advances
    the interior diffusion--recombination operator.
    """
    times = np.asarray(time_s, dtype=float)
    if times.ndim != 1 or times.size < 2 or times[0] != 0.0 or np.any(np.diff(times) <= 0.0):
        raise ValueError("time_s must increase strictly from zero")
    steady = minority_diffusion_fd(
        length_m, diffusivity_m2_s, lifetime_s, initial_junction_excess_m3, points=points
    )
    position = np.asarray(steady["position_m"])
    dx = position[1] - position[0]
    interior = points - 2
    laplacian = diags(
        (np.ones(interior - 1), -2.0 * np.ones(interior), np.ones(interior - 1)),
        (-1, 0, 1),
        format="csr",
    ) / dx**2
    operator = diffusivity_m2_s * laplacian - eye(interior, format="csr") / lifetime_s
    profile = np.zeros((times.size, points))
    # Retain the pre-switch interior profile while the switched contacts clamp
    # both boundary excess densities to zero at t=0+.
    profile[0, 1:-1] = np.asarray(steady["excess_m3"])[1:-1]
    for step, dt in enumerate(np.diff(times), start=1):
        profile[step, 1:-1] = spsolve(
            eye(interior, format="csr") - dt * operator,
            profile[step - 1, 1:-1],
        )
    stored_charge = Q * np.trapezoid(profile, position, axis=1)
    extraction_current = Q * diffusivity_m2_s * profile[:, 1] / dx
    return {
        "position_m": position,
        "time_s": times,
        "excess_m3": profile,
        "stored_charge_c_m2": stored_charge,
        "extraction_current_a_m2": extraction_current,
    }


def avalanche_integral(
    reverse_voltage_v: float,
    donor_m3: float,
    built_in_v: float,
    permittivity_f_m: float,
    prefactor_m1: float,
    characteristic_field_v_m: float,
    *,
    points: int = 1001,
) -> float:
    """Evaluate the equal-ionization-coefficient avalanche integral.

    Parameters
    ----------
    reverse_voltage_v : float
        Nonnegative reverse-bias magnitude in volts.
    donor_m3 : float
        Positive one-sided junction doping in ``m^-3``.
    built_in_v : float
        Positive built-in potential in volts.
    permittivity_f_m : float
        Positive permittivity in ``F m^-1``.
    prefactor_m1 : float
        Positive ionization prefactor ``A`` in ``m^-1``.
    characteristic_field_v_m : float
        Positive characteristic field ``B`` in ``V m^-1``.
    points : int, optional
        Number of quadrature samples across the depletion layer.

    Returns
    -------
    float
        Dimensionless integral ``integral alpha(x) dx``.

    Raises
    ------
    ValueError
        If parameters are invalid.
    """
    if reverse_voltage_v < 0.0 or min(donor_m3, built_in_v, permittivity_f_m, prefactor_m1, characteristic_field_v_m) <= 0.0:
        raise ValueError("bias must be nonnegative and parameters positive")
    total_drop = built_in_v + reverse_voltage_v
    width = np.sqrt(2.0 * permittivity_f_m * total_drop / (Q * donor_m3))
    maximum_field = Q * donor_m3 * width / permittivity_f_m
    position = np.linspace(0.0, width, points)
    field = maximum_field * (1.0 - position / width)
    alpha = np.zeros_like(field)
    positive = field > 0.0
    alpha[positive] = prefactor_m1 * np.exp(-characteristic_field_v_m / field[positive])
    return float(np.trapezoid(alpha, position))


def avalanche_breakdown_voltage(
    donor_m3: float,
    built_in_v: float,
    permittivity_f_m: float,
    prefactor_m1: float,
    characteristic_field_v_m: float,
    *,
    maximum_voltage_v: float = 5000.0,
) -> float:
    """Find reverse voltage at which the avalanche integral reaches unity.

    Parameters
    ----------
    donor_m3, built_in_v, permittivity_f_m : float
        Junction doping, built-in voltage, and permittivity in SI units.
    prefactor_m1, characteristic_field_v_m : float
        Positive impact-ionization parameters in SI units.
    maximum_voltage_v : float, optional
        Positive upper search limit in volts.

    Returns
    -------
    float
        Avalanche breakdown voltage in volts.

    Raises
    ------
    ValueError
        If the upper limit does not bracket breakdown.
    """
    function = lambda voltage: avalanche_integral(
        voltage,
        donor_m3,
        built_in_v,
        permittivity_f_m,
        prefactor_m1,
        characteristic_field_v_m,
    ) - 1.0
    if function(maximum_voltage_v) <= 0.0:
        raise ValueError("maximum_voltage_v does not bracket avalanche breakdown")
    return float(brentq(function, 0.0, maximum_voltage_v, xtol=1.0e-10))


def anderson_band_offsets(
    affinity_1_ev: float,
    gap_1_ev: float,
    affinity_2_ev: float,
    gap_2_ev: float,
) -> tuple[float, float]:
    """Calculate ideal heterojunction offsets with the electron-affinity rule.

    Parameters
    ----------
    affinity_1_ev : float
        Positive electron affinity of material 1 in electron-volts.
    gap_1_ev : float
        Positive band gap of material 1 in electron-volts.
    affinity_2_ev : float
        Positive electron affinity of material 2 in electron-volts.
    gap_2_ev : float
        Positive band gap of material 2 in electron-volts.

    Returns
    -------
    conduction_offset_ev : float
        Change ``E_C2-E_C1`` in electron-volts.
    valence_offset_ev : float
        Change ``E_V2-E_V1`` in electron-volts.

    Raises
    ------
    ValueError
        If an affinity or gap is nonpositive.
    """
    if min(affinity_1_ev, affinity_2_ev, gap_1_ev, gap_2_ev) <= 0.0:
        raise ValueError("affinities and gaps must be positive")
    conduction = affinity_1_ev - affinity_2_ev
    valence = affinity_1_ev + gap_1_ev - affinity_2_ev - gap_2_ev
    return conduction, valence


def equilibrium_junction_figure() -> None:
    """Write abrupt-junction depletion and Poisson--Boltzmann comparisons."""
    acceptor, donor, intrinsic = 5.0e22, 1.0e22, 1.0e16
    permittivity = 11.7 * EPSILON_0
    depletion = abrupt_depletion_profile(acceptor, donor, intrinsic, 300.0, permittivity)
    numerical = solve_poisson_boltzmann_junction(acceptor, donor, intrinsic, 300.0, permittivity)
    figure, axes = plt.subplots(1, 3, figsize=(11.0, 3.3), layout="constrained")
    x_dep = np.asarray(depletion["position_m"]) * 1.0e6
    axes[0].plot(x_dep, np.asarray(depletion["charge_c_m3"]) / 1.0e3)
    axes[0].set(xlabel="position (µm)", ylabel=r"charge ($10^3$ C m$^{-3}$)", title="Depletion charge")
    axes[1].plot(x_dep, np.asarray(depletion["field_v_m"]) / 1.0e5, label="depletion")
    axes[1].plot(np.asarray(numerical["position_m"]) * 1.0e6, np.asarray(numerical["field_v_m"]) / 1.0e5, "--", label="Poisson--Boltzmann")
    axes[1].set(xlabel="position (µm)", ylabel=r"field ($10^5$ V m$^{-1}$)", title="Junction field")
    axes[1].legend(fontsize=7)
    axes[2].plot(x_dep, np.asarray(depletion["potential_v"]), label="depletion")
    shifted = np.asarray(numerical["potential_v"]) - np.asarray(numerical["potential_v"])[0]
    axes[2].plot(np.asarray(numerical["position_m"]) * 1.0e6, shifted, "--", label="Poisson--Boltzmann")
    axes[2].set(xlabel="position (µm)", ylabel="potential (V)", title="Built-in potential")
    axes[2].legend(fontsize=7)
    for axis in axes:
        axis.grid(alpha=0.2)
    _save(figure, "junction-equilibrium-depletion.pdf")


def capacitance_figure() -> None:
    """Write depletion-capacitance and doping-profile validation curves."""
    acceptor, donor, intrinsic = 1.0e24, 2.0e21, 1.0e16
    permittivity = 11.7 * EPSILON_0
    reverse = np.linspace(0.0, 12.0, 301)
    capacitance = depletion_capacitance(reverse, acceptor, donor, intrinsic, 300.0, permittivity)
    extracted = extract_one_sided_doping(reverse, capacitance, permittivity)
    figure, axes = plt.subplots(1, 3, figsize=(11.0, 3.3), layout="constrained")
    # Convert F m^-2 to nF cm^-2: 1 F m^-2 = 10^5 nF cm^-2.
    axes[0].plot(reverse, capacitance * 1.0e5)
    axes[0].set(xlabel="reverse bias (V)", ylabel=r"$C_j$ (nF cm$^{-2}$)", title="Depletion capacitance")
    axes[1].plot(reverse, 1.0 / capacitance**2 / 1.0e14)
    axes[1].set(xlabel="reverse bias (V)", ylabel=r"$1/C_j^2$ ($10^{14}$ m$^4$ F$^{-2}$)", title="Abrupt-junction linearization")
    axes[2].plot(reverse, extracted / 1.0e6, label="extracted")
    axes[2].axhline(donor / 1.0e6, color="black", linestyle="--", label="input")
    axes[2].set(xlabel="reverse bias (V)", ylabel=r"doping (cm$^{-3}$)", title="C--V doping profile")
    axes[2].legend(fontsize=7)
    for axis in axes:
        axis.grid(alpha=0.2)
    _save(figure, "junction-capacitance-profiling.pdf")


def diode_current_figure() -> None:
    """Write ideal, recombination, and finite-base junction-current curves."""
    temperature = 300.0
    voltage = np.linspace(-0.2, 0.75, 500)
    ideal = ideal_diode_current_density(voltage, 1.0e-8, temperature)
    recombination = ideal_diode_current_density(voltage, 2.0e-5, temperature, ideality=2.0)
    length, diffusivity, lifetime = 120.0e-6, 1.0e-3, 2.0e-6
    diffusion_length = np.sqrt(diffusivity * lifetime)
    x = np.linspace(0.0, length, 201)
    numerical = minority_diffusion_fd(length, diffusivity, lifetime, 1.0, points=x.size)
    exact = np.sinh((length - x) / diffusion_length) / np.sinh(length / diffusion_length)
    figure, axes = plt.subplots(1, 3, figsize=(11.0, 3.3), layout="constrained")
    axes[0].semilogy(voltage, np.maximum(ideal, 1.0e-12), label="diffusion, η=1")
    axes[0].semilogy(voltage, np.maximum(recombination, 1.0e-12), label="depletion SRH, η=2")
    axes[0].set(xlabel="voltage (V)", ylabel=r"forward current density (A m$^{-2}$)", title="Ideality regimes")
    axes[0].legend(fontsize=7)
    axes[1].plot(x * 1.0e6, exact, label="analytic")
    axes[1].plot(x * 1.0e6, numerical["excess_m3"], "--", label="finite difference")
    axes[1].set(xlabel="neutral-region distance (µm)", ylabel="normalized excess", title="Minority diffusion")
    axes[1].legend(fontsize=7)
    bias = np.linspace(0.0, 0.7, 200)
    boundary = np.expm1(Q * bias / (K_B * temperature))
    flux = []
    for value in boundary:
        flux.append(minority_diffusion_fd(length, diffusivity, lifetime, value, points=101)["particle_flux_m2_s"])
    axes[2].semilogy(bias, np.maximum(Q * np.asarray(flux), 1.0e-20))
    axes[2].set(xlabel="forward bias (V)", ylabel=r"diffusion current density (A m$^{-2}$)", title="FD current extraction")
    for axis in axes:
        axis.grid(alpha=0.2, which="both")
    _save(figure, "junction-current-voltage.pdf")


def charge_storage_figure() -> None:
    """Write minority storage, extraction current, and diffusion-capacitance curves."""
    length, diffusivity, lifetime = 120.0e-6, 1.0e-3, 2.0e-6
    times = np.linspace(0.0, 12.0e-6, 121)
    transient = transient_minority_diffusion_fd(length, diffusivity, lifetime, 1.0e20, times, points=181)
    figure, axes = plt.subplots(1, 3, figsize=(11.0, 3.3), layout="constrained")
    for index in (0, 5, 20, 60, 120):
        axes[0].plot(
            np.asarray(transient["position_m"]) * 1.0e6,
            np.asarray(transient["excess_m3"])[index] / 1.0e20,
            label=fr"$t={times[index]*1e6:.1f}\,\mu$s",
        )
    axes[0].set(xlabel="position (µm)", ylabel=r"$\Delta p/10^{20}$ m$^{-3}$", title="Stored-charge removal")
    axes[0].legend(fontsize=6.5)
    charge = np.asarray(transient["stored_charge_c_m2"])
    axes[1].semilogy(times * 1.0e6, np.maximum(charge / charge[0], 1.0e-12))
    axes[1].set(xlabel="time (µs)", ylabel=r"$Q(t)/Q(0)$", title="Charge-storage decay")
    current = np.asarray(transient["extraction_current_a_m2"])
    axes[2].semilogy(times * 1.0e6, np.maximum(current, 1.0e-12))
    axes[2].set(xlabel="time (µs)", ylabel=r"extraction current (A m$^{-2}$)", title="Reverse-recovery current")
    for axis in axes:
        axis.grid(alpha=0.2, which="both")
    _save(figure, "junction-charge-storage.pdf")


def breakdown_heterojunction_figure() -> None:
    """Write avalanche-integral, breakdown-voltage, and band-offset curves."""
    permittivity = 11.7 * EPSILON_0
    built_in = 0.8
    voltage = np.linspace(0.0, 250.0, 300)
    dopings = (5.0e20, 1.0e21, 2.0e21)
    figure, axes = plt.subplots(1, 3, figsize=(11.0, 3.3), layout="constrained")
    for donor in dopings:
        integral = [avalanche_integral(v, donor, built_in, permittivity, 7.0e7, 1.23e8) for v in voltage]
        axes[0].semilogy(voltage, np.maximum(integral, 1.0e-12), label=fr"$N_D={donor/1e6:.0e}$ cm$^{{-3}}$")
    axes[0].axhline(1.0, color="black", linestyle="--")
    axes[0].set(xlabel="reverse bias (V)", ylabel=r"$\int\alpha\,dx$", title="Avalanche criterion")
    axes[0].legend(fontsize=6.5)
    doping_sweep = np.logspace(20.3, 22.5, 30)
    breakdown = [avalanche_breakdown_voltage(d, built_in, permittivity, 7.0e7, 1.23e8) for d in doping_sweep]
    axes[1].loglog(doping_sweep / 1.0e6, breakdown)
    axes[1].set(xlabel=r"light-side doping (cm$^{-3}$)", ylabel="breakdown voltage (V)", title="Doping controls breakdown")
    ec, ev = anderson_band_offsets(4.07, 1.42, 3.50, 2.16)
    x = np.array((-1.0, 0.0, 0.0, 1.0))
    conduction = np.array((0.0, 0.0, ec, ec))
    valence = np.array((-1.42, -1.42, ec - 2.16, ec - 2.16))
    axes[2].plot(x, conduction, label=r"$E_C$")
    axes[2].plot(x, valence, label=r"$E_V$")
    axes[2].axvline(0.0, color="0.5", linestyle=":")
    axes[2].set(xlabel="interface coordinate", ylabel="relative energy (eV)", title=fr"Heterojunction: $\Delta E_C={ec:.2f}$ eV")
    axes[2].legend()
    for axis in axes:
        axis.grid(alpha=0.2, which="both")
    _save(figure, "junction-breakdown-heterojunction.pdf")
