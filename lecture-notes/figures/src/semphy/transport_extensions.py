"""Reproducible figures for semiconductor transport beyond drift--diffusion.

The routines in this module support the Sze--Lee Chapter 2 sequence after the
microscopic Brownian, Drude--Sommerfeld, and effective-mass foundations.  They
cover mobility competition, Hall response, excess-carrier relaxation,
quasi-Fermi levels, barrier transport, space-charge-limited current, and
high-field effects.  Every numerical recipe uses explicit SI or declared
centimeter-based units and writes a timestamp-free PDF beneath ``figures/``.

The models are reduced teaching models.  They expose scaling laws and internal
consistency checks; they are not calibrated material or device simulators.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares, root
from scipy.sparse import eye, lil_matrix
from scipy.sparse.linalg import spsolve

FIGURE_DIR = Path(__file__).resolve().parents[2]
Q = 1.602_176_634e-19
K_B = 1.380_649e-23
HBAR = 1.054_571_817e-34
EPSILON_0 = 8.854_187_8128e-12
M_0 = 9.109_383_7139e-31

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


def save(figure: plt.Figure, name: str) -> None:
    """Write and close a deterministic publication figure.

    Parameters
    ----------
    figure : matplotlib.figure.Figure
        Fully configured Matplotlib figure to serialize.
    name : str
        Output filename beneath :data:`FIGURE_DIR`. The caller supplies the
        ``.pdf`` suffix.

    Returns
    -------
    None
        The function writes the PDF and closes ``figure``.

    Notes
    -----
    Creation and modification timestamps are omitted from the PDF metadata so
    identical inputs produce byte-stable publication artifacts.
    """
    # Strip time-varying metadata so rebuilding does not change the PDF bytes.
    metadata = {
        "Creator": "SEMPHY01 deterministic figure generator",
        "Producer": "Matplotlib",
        "CreationDate": None,
        "ModDate": None,
    }
    figure.savefig(FIGURE_DIR / name, metadata=metadata)
    plt.close(figure)


def matthiessen_mobility(
    lattice_mobility: np.ndarray | float,
    impurity_mobility: np.ndarray | float,
) -> np.ndarray:
    """Combine lattice and impurity mobility limits with Matthiessen's rule.

    Parameters
    ----------
    lattice_mobility : array_like
        Positive lattice-limited mobility. Values may be in SI units or
        ``cm^2 V^-1 s^-1``.
    impurity_mobility : array_like
        Positive ionized-impurity-limited mobility in the same unit as
        ``lattice_mobility``. Inputs follow NumPy broadcasting rules.

    Returns
    -------
    numpy.ndarray
        Total mobility in the common input unit.

    Raises
    ------
    ValueError
        If either mobility contains a nonpositive value.

    Notes
    -----
    Independent momentum-relaxation rates give
    ``1 / mu = 1 / mu_lattice + 1 / mu_impurity``. The result is therefore
    smaller than either individual mobility limit.
    """
    lattice = np.asarray(lattice_mobility, dtype=float)
    impurity = np.asarray(impurity_mobility, dtype=float)
    if np.any(lattice <= 0.0) or np.any(impurity <= 0.0):
        raise ValueError("mobility limits must be positive")
    # Independent scattering rates add, so inverse mobilities add as well.
    return 1.0 / (1.0 / lattice + 1.0 / impurity)


def hall_coefficient(density_m3: np.ndarray | float, charge_sign: int) -> np.ndarray:
    """Calculate the ideal one-carrier Hall coefficient.

    Parameters
    ----------
    density_m3 : array_like
        Positive mobile-carrier density in ``m^-3``.
    charge_sign : {-1, 1}
        Carrier-charge sign: ``-1`` for electrons and ``+1`` for holes.

    Returns
    -------
    numpy.ndarray
        Hall coefficient in ``m^3 C^-1`` with the broadcast shape of
        ``density_m3``.

    Raises
    ------
    ValueError
        If a density is nonpositive or ``charge_sign`` is not ``-1`` or ``1``.

    Notes
    -----
    The result is ``R_H = charge_sign / (q N)``. It assumes one isotropic
    carrier population and a Hall scattering factor of unity.
    """
    density = np.asarray(density_m3, dtype=float)
    if np.any(density <= 0.0):
        raise ValueError("carrier density must be positive")
    if charge_sign not in (-1, 1):
        raise ValueError("charge_sign must be -1 for electrons or +1 for holes")
    # q is stored as a positive magnitude; charge_sign carries the Hall polarity.
    return charge_sign / (Q * density)


def two_carrier_hall_coefficient(
    electron_density_m3: np.ndarray | float,
    hole_density_m3: np.ndarray | float,
    electron_mobility: np.ndarray | float,
    hole_mobility: np.ndarray | float,
) -> np.ndarray:
    """Calculate the weak-field Hall coefficient for electrons and holes.

    Parameters
    ----------
    electron_density_m3 : array_like
        Nonnegative electron density ``n`` in ``m^-3``.
    hole_density_m3 : array_like
        Nonnegative hole density ``p`` in ``m^-3``.
    electron_mobility : array_like
        Positive electron mobility ``mu_n``.
    hole_mobility : array_like
        Positive hole mobility ``mu_p`` in the same unit as
        ``electron_mobility``. All inputs follow NumPy broadcasting rules.

    Returns
    -------
    numpy.ndarray
        Weak-field Hall coefficient in ``m^3 C^-1``.

    Raises
    ------
    ValueError
        If a density is negative, a mobility is nonpositive, or both carrier
        populations provide zero conductivity.

    Notes
    -----
    The implemented expression is
    ``(p*mu_p**2 - n*mu_n**2) / (q*(p*mu_p + n*mu_n)**2)``. Mobilities may use
    any common unit because their conversion factor cancels. The squared
    mobility weighting allows a minority population to control the Hall sign.
    """
    electron = np.asarray(electron_density_m3, dtype=float)
    hole = np.asarray(hole_density_m3, dtype=float)
    mu_e = np.asarray(electron_mobility, dtype=float)
    mu_h = np.asarray(hole_mobility, dtype=float)
    if np.any(electron < 0.0) or np.any(hole < 0.0):
        raise ValueError("carrier densities must be nonnegative")
    if np.any(mu_e <= 0.0) or np.any(mu_h <= 0.0):
        raise ValueError("mobilities must be positive")
    # The denominator is proportional to the square of zero-field conductivity.
    conductivity_weight = electron * mu_e + hole * mu_h
    if np.any(conductivity_weight <= 0.0):
        raise ValueError("at least one carrier population must be positive")
    # Opposite charge signs enter the numerator with mobility-squared weights.
    return (hole * mu_h**2 - electron * mu_e**2) / (
        Q * conductivity_weight**2
    )


def excess_decay(
    time_s: np.ndarray | float,
    initial_excess_m3: float,
    lifetime_s: float,
) -> np.ndarray:
    """Calculate low-injection excess-carrier decay after generation stops.

    Parameters
    ----------
    time_s : array_like
        Nonnegative time after the excitation is removed, in seconds.
    initial_excess_m3 : float
        Nonnegative initial excess-carrier density in ``m^-3``.
    lifetime_s : float
        Positive minority-carrier lifetime in seconds.

    Returns
    -------
    numpy.ndarray
        Excess density ``Delta n(0) exp(-t/tau)`` in ``m^-3``.

    Raises
    ------
    ValueError
        If time or initial excess is negative, or if the lifetime is not
        positive.
    """
    time = np.asarray(time_s, dtype=float)
    if np.any(time < 0.0) or initial_excess_m3 < 0.0 or lifetime_s <= 0.0:
        raise ValueError("time and excess must be nonnegative and lifetime positive")
    return initial_excess_m3 * np.exp(-time / lifetime_s)


def quasi_fermi_splitting(
    electron_m3: np.ndarray | float,
    hole_m3: np.ndarray | float,
    intrinsic_m3: float,
    temperature_k: float,
) -> np.ndarray:
    """Calculate the nondegenerate quasi-Fermi-level splitting.

    Parameters
    ----------
    electron_m3 : array_like
        Positive electron density in ``m^-3``.
    hole_m3 : array_like
        Positive hole density in ``m^-3``.
    intrinsic_m3 : float
        Positive intrinsic carrier density in ``m^-3``.
    temperature_k : float
        Positive absolute temperature in kelvin.

    Returns
    -------
    numpy.ndarray
        Splitting ``E_Fn - E_Fp`` in electron-volts.

    Raises
    ------
    ValueError
        If any carrier density or the temperature is nonpositive.

    Notes
    -----
    The Boltzmann relation is
    ``E_Fn - E_Fp = k_B T ln(n p / n_i**2)``. Dividing ``k_B T`` in joules by
    the elementary charge gives the numerical energy in electron-volts.
    """
    electron = np.asarray(electron_m3, dtype=float)
    hole = np.asarray(hole_m3, dtype=float)
    if np.any(electron <= 0.0) or np.any(hole <= 0.0):
        raise ValueError("carrier densities must be positive")
    if intrinsic_m3 <= 0.0 or temperature_k <= 0.0:
        raise ValueError("intrinsic density and temperature must be positive")
    # k_B T / q is numerically both a voltage and an energy in electron-volts.
    thermal_voltage = K_B * temperature_k / Q
    return thermal_voltage * np.log(electron * hole / intrinsic_m3**2)


def rectangular_barrier_transmission(
    energy_j: np.ndarray | float,
    barrier_j: float,
    width_m: float,
    effective_mass_kg: float,
) -> np.ndarray:
    """Calculate transmission through a one-dimensional rectangular barrier.

    Parameters
    ----------
    energy_j : array_like
        Incident particle energy in joules, restricted to ``0 < E < U0``.
    barrier_j : float
        Positive barrier energy ``U0`` in joules.
    width_m : float
        Positive barrier width in meters.
    effective_mass_kg : float
        Positive carrier effective mass in kilograms.

    Returns
    -------
    numpy.ndarray
        Dimensionless transmission probability between zero and one.

    Raises
    ------
    ValueError
        If an energy lies outside the sub-barrier interval or if a model
        parameter is nonpositive.

    Notes
    -----
    This is the exact constant-mass result for a rectangular barrier. The
    argument of ``sinh`` is clipped at 350 to prevent floating-point overflow
    in the exponentially opaque limit without changing the resolved result.
    """
    energy = np.asarray(energy_j, dtype=float)
    if np.any(energy <= 0.0) or np.any(energy >= barrier_j):
        raise ValueError("energy must satisfy 0 < energy < barrier")
    if barrier_j <= 0.0 or width_m <= 0.0 or effective_mass_kg <= 0.0:
        raise ValueError("barrier, width, and effective mass must be positive")
    # Evanescent decay constant inside the classically forbidden barrier.
    kappa = np.sqrt(2.0 * effective_mass_kg * (barrier_j - energy)) / HBAR
    denominator = 1.0 + (
        barrier_j**2
        # Clipping avoids overflow only after transmission is numerically opaque.
        * np.sinh(np.minimum(kappa * width_m, 350.0)) ** 2
        / (4.0 * energy * (barrier_j - energy))
    )
    return 1.0 / denominator


def fd_rectangular_barrier_transmission(
    energy_j: np.ndarray | float,
    barrier_j: float,
    width_m: float,
    effective_mass_kg: float,
    grid_points: int,
) -> np.ndarray:
    """Calculate rectangular-barrier transmission on a finite-difference grid.

    Parameters
    ----------
    energy_j : array_like
        Positive incident energy in joules below ``barrier_j``.
    barrier_j : float
        Positive rectangular-barrier energy in joules.
    width_m : float
        Positive barrier width in meters.
    effective_mass_kg : float
        Positive constant effective mass in kilograms.
    grid_points : int
        Number of cell-centered finite-difference sites across the barrier.

    Returns
    -------
    numpy.ndarray
        Dimensionless transmission probability with the shape of ``energy_j``.

    Raises
    ------
    ValueError
        If energies are outside ``0 < E < barrier_j``, a physical parameter is
        nonpositive, or fewer than four grid sites are requested.

    Notes
    -----
    The second derivative produces a nearest-neighbor Hamiltonian with hopping
    ``t0 = hbar**2/(2*m*dx**2)``. Semi-infinite zero-potential leads enter
    through analytic discrete self-energies. The Caroli expression reduces to
    ``T = Gamma_L*Gamma_R*abs(G_1N)**2`` for the one-dimensional chain.
    """
    energy = np.asarray(energy_j, dtype=float)
    if np.any(energy <= 0.0) or np.any(energy >= barrier_j):
        raise ValueError("energy must satisfy 0 < energy < barrier")
    if barrier_j <= 0.0 or width_m <= 0.0 or effective_mass_kg <= 0.0:
        raise ValueError("barrier, width, and effective mass must be positive")
    if type(grid_points) is not int or grid_points < 4:
        raise ValueError("grid_points must be an integer of at least four")

    dx = width_m / grid_points
    hopping = HBAR**2 / (2.0 * effective_mass_kg * dx**2)
    transmission = np.empty_like(energy)
    flat_energy = energy.ravel()
    flat_transmission = transmission.ravel()

    # The cell-centered device contains the barrier; both leads have V=0.
    for index, value in enumerate(flat_energy):
        cosine = 1.0 - value / (2.0 * hopping)
        if not -1.0 < cosine < 1.0:
            raise ValueError("energy lies outside the propagating lead band")
        wave_number_dx = np.arccos(cosine)
        self_energy = -hopping * np.exp(1j * wave_number_dx)
        broadening = -2.0 * np.imag(self_energy)

        system = lil_matrix((grid_points, grid_points), dtype=complex)
        system.setdiag(value - (2.0 * hopping + barrier_j))
        system.setdiag(np.full(grid_points - 1, hopping), k=1)
        system.setdiag(np.full(grid_points - 1, hopping), k=-1)
        # Subtracting the retarded lead self-energy opens the finite device.
        system[0, 0] -= self_energy
        system[-1, -1] -= self_energy
        source = np.zeros(grid_points, dtype=complex)
        source[-1] = 1.0
        endpoint_green = spsolve(system.tocsr(), source)[0]
        flat_transmission[index] = broadening**2 * abs(endpoint_green) ** 2
    return transmission


def poisson_potential_fd(
    position_m: np.ndarray,
    charge_density_c_m3: np.ndarray | float,
    permittivity_f_m: float,
    left_potential_v: float,
    right_potential_v: float,
) -> np.ndarray:
    """Solve one-dimensional Poisson's equation with Dirichlet boundaries.

    Parameters
    ----------
    position_m : numpy.ndarray
        Uniform grid in meters including both boundary points.
    charge_density_c_m3 : array_like
        Charge density ``rho(x)`` in ``C m^-3``.
    permittivity_f_m : float
        Positive permittivity in ``F m^-1``.
    left_potential_v, right_potential_v : float
        Electrostatic potentials at the two boundaries in volts.

    Returns
    -------
    numpy.ndarray
        Electrostatic potential in volts on ``position_m``.

    Raises
    ------
    ValueError
        If the grid is nonuniform, dimensions disagree, or permittivity is not
        positive.

    Notes
    -----
    Interior rows discretize ``d2phi/dx2 = -rho/epsilon`` with centered second
    differences. The boundary rows impose the two contact potentials exactly.
    """
    position = np.asarray(position_m, dtype=float)
    if position.ndim != 1 or position.size < 3:
        raise ValueError("position_m must contain at least three points")
    spacing = np.diff(position)
    if not np.allclose(spacing, spacing[0], rtol=1.0e-10, atol=0.0):
        raise ValueError("position_m must be uniformly spaced")
    charge = np.broadcast_to(
        np.asarray(charge_density_c_m3, dtype=float), position.shape
    ).copy()
    if permittivity_f_m <= 0.0:
        raise ValueError("permittivity_f_m must be positive")

    dx = spacing[0]
    interior_count = position.size - 2
    matrix = lil_matrix((interior_count, interior_count), dtype=float)
    matrix.setdiag(np.full(interior_count, -2.0))
    if interior_count > 1:
        matrix.setdiag(np.ones(interior_count - 1), k=1)
        matrix.setdiag(np.ones(interior_count - 1), k=-1)
    right_hand_side = -dx**2 * charge[1:-1] / permittivity_f_m
    right_hand_side[0] -= left_potential_v
    right_hand_side[-1] -= right_potential_v
    potential = np.empty_like(position)
    potential[0] = left_potential_v
    potential[-1] = right_potential_v
    potential[1:-1] = spsolve(matrix.tocsr(), right_hand_side)
    return potential


def space_charge_field_fd(
    voltage_v: float,
    length_m: float,
    permittivity_f_m: float,
    grid_points: int,
    *,
    mobility_m2_v_s: float | None = None,
    saturation_velocity_m_s: float | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Solve a one-dimensional space-charge-limited field by finite differences.

    Parameters
    ----------
    voltage_v : float
        Positive applied voltage in volts.
    length_m : float
        Positive layer length in meters.
    permittivity_f_m : float
        Positive semiconductor permittivity in ``F m^-1``.
    grid_points : int
        Number of field nodes including both contacts.
    mobility_m2_v_s : float, optional
        Constant positive mobility in ``m^2 V^-1 s^-1``.
    saturation_velocity_m_s : float, optional
        Constant positive drift speed in ``m s^-1``.

    Returns
    -------
    position_m : numpy.ndarray
        Uniform field grid in meters.
    field_v_m : numpy.ndarray
        Nonnegative electric field in ``V m^-1``.
    current_density_a_m2 : float
        Space-charge-limited current density in ``A m^-2``.

    Raises
    ------
    ValueError
        If parameters are invalid or if both/neither transport laws are given.
    RuntimeError
        If the nonlinear finite-difference solve does not converge.

    Notes
    -----
    Exactly one transport law is selected. Constant mobility uses
    ``v(E)=mu*E`` and constant velocity uses ``v(E)=v_sat``. Cell residuals
    discretize ``dE/dx=J/(epsilon*v(E))``. The injecting-contact condition is
    ``E(0)=0`` and a trapezoidal voltage constraint determines the unknown
    current density.
    """
    if voltage_v <= 0.0 or length_m <= 0.0 or permittivity_f_m <= 0.0:
        raise ValueError("voltage, length, and permittivity must be positive")
    if type(grid_points) is not int or grid_points < 4:
        raise ValueError("grid_points must be an integer of at least four")
    selected = (mobility_m2_v_s is not None) + (saturation_velocity_m_s is not None)
    if selected != 1:
        raise ValueError("provide exactly one of mobility or saturation velocity")
    if mobility_m2_v_s is not None and mobility_m2_v_s <= 0.0:
        raise ValueError("mobility must be positive")
    if saturation_velocity_m_s is not None and saturation_velocity_m_s <= 0.0:
        raise ValueError("saturation velocity must be positive")

    position = np.linspace(0.0, length_m, grid_points)
    dx = position[1] - position[0]
    if mobility_m2_v_s is not None:
        analytic_current = (
            9.0 * permittivity_f_m * mobility_m2_v_s * voltage_v**2
            / (8.0 * length_m**3)
        )
        initial_field = np.sqrt(
            2.0 * analytic_current * position
            / (permittivity_f_m * mobility_m2_v_s)
        )
    else:
        analytic_current = (
            2.0 * permittivity_f_m * saturation_velocity_m_s * voltage_v
            / length_m**2
        )
        initial_field = 2.0 * voltage_v * position / length_m**2
    field_scale = voltage_v / length_m
    initial = np.concatenate((initial_field / field_scale, (analytic_current,)))

    def normalized_residual(unknown: np.ndarray) -> np.ndarray:
        field = unknown[:-1] * field_scale
        current = unknown[-1]
        midpoint_field = 0.5 * (field[:-1] + field[1:])
        if mobility_m2_v_s is not None:
            speed = mobility_m2_v_s * midpoint_field
        else:
            speed = np.full_like(midpoint_field, saturation_velocity_m_s)
        residual = np.empty(grid_points + 1)
        residual[0] = field[0] / field_scale
        residual[1:-1] = (
            field[1:]
            - field[:-1]
            - dx * current / (permittivity_f_m * speed)
        ) / field_scale
        residual[-1] = (np.trapezoid(field, position) - voltage_v) / voltage_v
        return residual

    lower = np.concatenate((np.zeros(grid_points), (0.0,)))
    upper = np.full(grid_points + 1, np.inf)
    result = least_squares(
        normalized_residual,
        initial,
        bounds=(lower, upper),
        xtol=1.0e-12,
        ftol=1.0e-12,
        gtol=1.0e-12,
        max_nfev=4000,
    )
    if not result.success or np.max(np.abs(normalized_residual(result.x))) > 1.0e-7:
        raise RuntimeError("space-charge finite-difference solve did not converge")
    return position, result.x[:-1] * field_scale, float(result.x[-1])


def mott_gurney_current(
    voltage_v: np.ndarray | float,
    permittivity_f_m: float,
    mobility_m2_v_s: float,
    length_m: float,
) -> np.ndarray:
    """Calculate trap-free constant-mobility space-charge-limited current.

    Parameters
    ----------
    voltage_v : array_like
        Nonnegative applied voltage in volts.
    permittivity_f_m : float
        Positive material permittivity in ``F m^-1``.
    mobility_m2_v_s : float
        Positive carrier mobility in ``m^2 V^-1 s^-1``.
    length_m : float
        Positive transport length in meters.

    Returns
    -------
    numpy.ndarray
        Mott--Gurney current density in ``A m^-2``.

    Raises
    ------
    ValueError
        If voltage is negative or a material parameter is nonpositive.

    Notes
    -----
    The model evaluates ``J = 9 epsilon mu V**2 / (8 L**3)`` for an ohmic
    injecting contact, constant mobility, and no traps.
    """
    voltage = np.asarray(voltage_v, dtype=float)
    if np.any(voltage < 0.0):
        raise ValueError("voltage must be nonnegative")
    if permittivity_f_m <= 0.0 or mobility_m2_v_s <= 0.0 or length_m <= 0.0:
        raise ValueError("permittivity, mobility, and length must be positive")
    return 9.0 * permittivity_f_m * mobility_m2_v_s * voltage**2 / (8.0 * length_m**3)


def saturated_space_charge_current(
    voltage_v: np.ndarray | float,
    permittivity_f_m: float,
    saturation_velocity_m_s: float,
    length_m: float,
) -> np.ndarray:
    """Calculate velocity-saturated space-charge-limited current density.

    Parameters
    ----------
    voltage_v : array_like
        Nonnegative applied voltage in volts.
    permittivity_f_m : float
        Positive material permittivity in ``F m^-1``.
    saturation_velocity_m_s : float
        Positive saturation velocity in ``m s^-1``.
    length_m : float
        Positive transport length in meters.

    Returns
    -------
    numpy.ndarray
        Velocity-saturated current density in ``A m^-2``.

    Raises
    ------
    ValueError
        If voltage is negative or a material parameter is nonpositive.

    Notes
    -----
    Once the drift speed is fixed at ``v_sat``, the reduced space-charge model
    becomes ``J = 2 epsilon v_sat V / L**2`` rather than the quadratic
    Mott--Gurney law.
    """
    voltage = np.asarray(voltage_v, dtype=float)
    if np.any(voltage < 0.0):
        raise ValueError("voltage must be nonnegative")
    if permittivity_f_m <= 0.0 or saturation_velocity_m_s <= 0.0 or length_m <= 0.0:
        raise ValueError("permittivity, velocity, and length must be positive")
    return 2.0 * permittivity_f_m * saturation_velocity_m_s * voltage / length_m**2


def saturation_velocity_model(
    field_v_m: np.ndarray | float,
    mobility_m2_v_s: float,
    saturation_velocity_m_s: float,
    exponent: float,
) -> np.ndarray:
    """Interpolate between low-field drift and velocity saturation.

    Parameters
    ----------
    field_v_m : array_like
        Nonnegative electric-field magnitude in ``V m^-1``.
    mobility_m2_v_s : float
        Positive low-field mobility in ``m^2 V^-1 s^-1``.
    saturation_velocity_m_s : float
        Positive limiting drift speed in ``m s^-1``.
    exponent : float
        Positive empirical sharpness parameter.

    Returns
    -------
    numpy.ndarray
        Drift-speed magnitude in ``m s^-1``.

    Raises
    ------
    ValueError
        If field is negative or a model parameter is nonpositive.

    Notes
    -----
    The interpolation approaches ``mu E`` when ``mu E << v_sat`` and
    ``v_sat`` when ``mu E >> v_sat``. The exponent controls the crossover but
    does not alter either asymptote.
    """
    field = np.asarray(field_v_m, dtype=float)
    if np.any(field < 0.0):
        raise ValueError("field must be nonnegative")
    if mobility_m2_v_s <= 0.0 or saturation_velocity_m_s <= 0.0 or exponent <= 0.0:
        raise ValueError("mobility, saturation velocity, and exponent must be positive")
    # reduced=1 marks the crossover mu*E approximately equal to v_sat.
    reduced = mobility_m2_v_s * field / saturation_velocity_m_s
    return mobility_m2_v_s * field / (1.0 + reduced**exponent) ** (1.0 / exponent)


def ionization_rate(
    field_v_m: np.ndarray | float,
    prefactor_m1: float,
    characteristic_field_v_m: float,
) -> np.ndarray:
    """Calculate an empirical impact-ionization coefficient.

    Parameters
    ----------
    field_v_m : array_like
        Positive electric-field magnitude in ``V m^-1``.
    prefactor_m1 : float
        Positive ionization prefactor ``A`` in ``m^-1``.
    characteristic_field_v_m : float
        Positive characteristic field ``B`` in ``V m^-1``.

    Returns
    -------
    numpy.ndarray
        Impact-ionization coefficient ``alpha`` in ``m^-1``.

    Raises
    ------
    ValueError
        If the field or either empirical parameter is nonpositive.

    Notes
    -----
    The model evaluates ``alpha(E) = A exp(-B/E)`` and captures the rapid
    field activation without modeling the microscopic hot-carrier distribution.
    """
    field = np.asarray(field_v_m, dtype=float)
    if np.any(field <= 0.0):
        raise ValueError("field must be positive")
    if prefactor_m1 <= 0.0 or characteristic_field_v_m <= 0.0:
        raise ValueError("ionization parameters must be positive")
    return prefactor_m1 * np.exp(-characteristic_field_v_m / field)


def srh_recombination_rate(
    electron_m3: np.ndarray | float,
    hole_m3: np.ndarray | float,
    intrinsic_m3: float,
    electron_lifetime_s: float,
    hole_lifetime_s: float,
    electron_auxiliary_m3: float,
    hole_auxiliary_m3: float,
) -> np.ndarray:
    """Calculate the Shockley--Read--Hall net recombination rate.

    Parameters
    ----------
    electron_m3 : array_like
        Positive electron density ``n`` in ``m^-3``.
    hole_m3 : array_like
        Positive hole density ``p`` in ``m^-3``.
    intrinsic_m3 : float
        Positive intrinsic density ``n_i`` in ``m^-3``.
    electron_lifetime_s : float
        Positive electron capture lifetime ``tau_n`` in seconds.
    hole_lifetime_s : float
        Positive hole capture lifetime ``tau_p`` in seconds.
    electron_auxiliary_m3 : float
        Positive SRH auxiliary density ``n_B`` in ``m^-3``.
    hole_auxiliary_m3 : float
        Positive SRH auxiliary density ``p_B`` in ``m^-3``.

    Returns
    -------
    numpy.ndarray
        Net SRH rate in ``m^-3 s^-1``. Positive values denote net
        recombination and negative values denote net generation.

    Raises
    ------
    ValueError
        If a density or lifetime is nonpositive.

    Notes
    -----
    The auxiliary densities include the center-degeneracy convention and must
    satisfy ``n_B*p_B = n_i**2`` for equilibrium detailed balance.
    """
    electron = np.asarray(electron_m3, dtype=float)
    hole = np.asarray(hole_m3, dtype=float)
    positive_parameters = (
        intrinsic_m3,
        electron_lifetime_s,
        hole_lifetime_s,
        electron_auxiliary_m3,
        hole_auxiliary_m3,
    )
    if np.any(electron <= 0.0) or np.any(hole <= 0.0):
        raise ValueError("carrier densities must be positive")
    if any(value <= 0.0 for value in positive_parameters):
        raise ValueError("intrinsic, lifetime, and auxiliary values must be positive")

    denominator = (
        hole_lifetime_s * (electron + electron_auxiliary_m3)
        + electron_lifetime_s * (hole + hole_auxiliary_m3)
    )
    return (electron * hole - intrinsic_m3**2) / denominator


def transient_continuity_fd(
    position_m: np.ndarray,
    time_s: np.ndarray,
    initial_excess_m3: np.ndarray,
    diffusivity_m2_s: float,
    drift_velocity_m_s: float,
    lifetime_s: float,
    generation_m3_s: np.ndarray | float = 0.0,
    *,
    theta: float = 0.5,
) -> np.ndarray:
    """Solve a periodic one-dimensional continuity equation by finite differences.

    Parameters
    ----------
    position_m : numpy.ndarray
        Uniform, endpoint-excluded periodic grid in meters.
    time_s : numpy.ndarray
        Strictly increasing output times in seconds. The initial condition is
        assigned at ``time_s[0]``.
    initial_excess_m3 : numpy.ndarray
        Initial excess-carrier density on ``position_m`` in ``m^-3``.
    diffusivity_m2_s : float
        Positive diffusivity ``D`` in ``m^2 s^-1``.
    drift_velocity_m_s : float
        Signed drift velocity ``v_d`` in ``m s^-1``.
    lifetime_s : float
        Positive linear recombination lifetime in seconds.
    generation_m3_s : array_like, optional
        Time-independent generation rate on the grid in ``m^-3 s^-1``.
    theta : float, optional
        Implicitness parameter. ``0.5`` gives Crank--Nicolson and ``1`` gives
        backward Euler.

    Returns
    -------
    numpy.ndarray
        Excess density with shape ``(time_s.size, position_m.size)``.

    Raises
    ------
    ValueError
        If a grid is nonuniform, dimensions disagree, times do not increase,
        a physical coefficient is invalid, or ``theta`` lies outside
        ``[0.5, 1]``.

    Notes
    -----
    Centered differences represent diffusion and second-order upwind
    differences represent drift. The theta method advances
    ``du/dt = D*d2u/dx2 - v_d*du/dx - u/tau + G``. Periodic boundaries isolate
    discretization error when comparing a moving packet with the analytic
    Haynes--Shockley solution.
    """
    position = np.asarray(position_m, dtype=float)
    times = np.asarray(time_s, dtype=float)
    initial = np.asarray(initial_excess_m3, dtype=float)
    if position.ndim != 1 or position.size < 4:
        raise ValueError("position_m must contain at least four grid points")
    if times.ndim != 1 or times.size < 2 or np.any(np.diff(times) <= 0.0):
        raise ValueError("time_s must be a strictly increasing one-dimensional grid")
    if initial.shape != position.shape or np.any(initial < 0.0):
        raise ValueError("initial_excess_m3 must be nonnegative and match position_m")
    spacing = np.diff(position)
    if not np.allclose(spacing, spacing[0], rtol=1.0e-10, atol=0.0):
        raise ValueError("position_m must be uniformly spaced")
    if diffusivity_m2_s <= 0.0 or lifetime_s <= 0.0:
        raise ValueError("diffusivity and lifetime must be positive")
    if not 0.5 <= theta <= 1.0:
        raise ValueError("theta must lie between 0.5 and 1")

    generation = np.broadcast_to(
        np.asarray(generation_m3_s, dtype=float), position.shape
    ).copy()
    if np.any(generation < 0.0):
        raise ValueError("generation_m3_s must be nonnegative")

    dx = spacing[0]
    point_count = position.size
    operator = lil_matrix((point_count, point_count), dtype=float)
    diffusion_weight = diffusivity_m2_s / dx**2
    drift_weight = abs(drift_velocity_m_s) / dx

    # Periodic indexing keeps a translating packet away from artificial contacts.
    for index in range(point_count):
        operator[index, index] = (
            -2.0 * diffusion_weight - 1.5 * drift_weight - 1.0 / lifetime_s
        )
        operator[index, (index - 1) % point_count] += diffusion_weight
        operator[index, (index + 1) % point_count] += diffusion_weight
        # A second-order upwind stencil suppresses the artificial packet
        # broadening produced by the first-order donor-cell formula.
        if drift_velocity_m_s >= 0.0:
            operator[index, (index - 1) % point_count] += 2.0 * drift_weight
            operator[index, (index - 2) % point_count] -= 0.5 * drift_weight
        else:
            operator[index, (index + 1) % point_count] += 2.0 * drift_weight
            operator[index, (index + 2) % point_count] -= 0.5 * drift_weight
    operator = operator.tocsr()
    identity = eye(point_count, format="csr")

    solution = np.empty((times.size, point_count), dtype=float)
    solution[0] = initial
    for step, dt in enumerate(np.diff(times), start=1):
        left = identity - theta * dt * operator
        right = (identity + (1.0 - theta) * dt * operator) @ solution[step - 1]
        right += dt * generation
        solution[step] = spsolve(left, right)
    return solution


def steady_srh_profile_fd(
    position_m: np.ndarray,
    generation_m3_s: np.ndarray | float,
    diffusivity_m2_s: float,
    equilibrium_electron_m3: float,
    equilibrium_hole_m3: float,
    intrinsic_m3: float,
    electron_lifetime_s: float,
    hole_lifetime_s: float,
    electron_auxiliary_m3: float,
    hole_auxiliary_m3: float,
    surface_velocity_m_s: float,
    *,
    right_excess_m3: float = 0.0,
) -> np.ndarray:
    """Solve a nonlinear steady minority-carrier profile by finite differences.

    Parameters
    ----------
    position_m : numpy.ndarray
        Uniform one-dimensional grid from the recombining surface to the right
        contact, in meters.
    generation_m3_s : array_like
        Nonnegative pair-generation rate in ``m^-3 s^-1``.
    diffusivity_m2_s : float
        Positive minority-carrier diffusivity in ``m^2 s^-1``.
    equilibrium_electron_m3, equilibrium_hole_m3 : float
        Positive equilibrium carrier densities in ``m^-3``.
    intrinsic_m3 : float
        Positive intrinsic density in ``m^-3``.
    electron_lifetime_s, hole_lifetime_s : float
        Positive SRH capture lifetimes in seconds.
    electron_auxiliary_m3, hole_auxiliary_m3 : float
        Positive SRH auxiliary densities in ``m^-3``.
    surface_velocity_m_s : float
        Nonnegative recombination velocity at the left surface in ``m s^-1``.
    right_excess_m3 : float, optional
        Prescribed nonnegative excess density at the right boundary.

    Returns
    -------
    numpy.ndarray
        Nonnegative quasineutral excess density ``Delta n = Delta p`` in
        ``m^-3`` on the grid.

    Raises
    ------
    ValueError
        If the grid, dimensions, or physical parameters are invalid.
    RuntimeError
        If the nonlinear root solve does not converge to a nonnegative profile.

    Notes
    -----
    Interior nodes satisfy ``D*d2Delta/dx2 + G - U_SRH = 0``. The left Robin
    condition is ``D*dDelta/dx = S*Delta`` and the right density is prescribed.
    The SRH rate is evaluated with ``n=n0+Delta`` and ``p=p0+Delta`` rather than
    replaced by a constant lifetime.
    """
    position = np.asarray(position_m, dtype=float)
    if position.ndim != 1 or position.size < 4:
        raise ValueError("position_m must contain at least four grid points")
    spacing = np.diff(position)
    if not np.allclose(spacing, spacing[0], rtol=1.0e-10, atol=0.0):
        raise ValueError("position_m must be uniformly spaced")
    generation = np.broadcast_to(
        np.asarray(generation_m3_s, dtype=float), position.shape
    ).copy()
    positive_parameters = (
        diffusivity_m2_s,
        equilibrium_electron_m3,
        equilibrium_hole_m3,
        intrinsic_m3,
        electron_lifetime_s,
        hole_lifetime_s,
        electron_auxiliary_m3,
        hole_auxiliary_m3,
    )
    if np.any(generation < 0.0):
        raise ValueError("generation_m3_s must be nonnegative")
    if any(value <= 0.0 for value in positive_parameters):
        raise ValueError("diffusivity, densities, and lifetimes must be positive")
    if surface_velocity_m_s < 0.0 or right_excess_m3 < 0.0:
        raise ValueError("surface velocity and boundary excess must be nonnegative")

    dx = spacing[0]
    source_scale = np.max(generation) * max(electron_lifetime_s, hole_lifetime_s)
    concentration_scale = max(source_scale, right_excess_m3, 1.0)
    initial_level = max(source_scale, right_excess_m3)
    initial = np.linspace(initial_level, right_excess_m3, position.size)

    def normalized_residual(excess: np.ndarray) -> np.ndarray:
        electron = equilibrium_electron_m3 + excess
        hole = equilibrium_hole_m3 + excess
        recombination = srh_recombination_rate(
            electron,
            hole,
            intrinsic_m3,
            electron_lifetime_s,
            hole_lifetime_s,
            electron_auxiliary_m3,
            hole_auxiliary_m3,
        )
        residual = np.empty_like(excess)
        # Multiplication by dx**2/D puts interior equations in density units.
        residual[1:-1] = (
            excess[:-2]
            - 2.0 * excess[1:-1]
            + excess[2:]
            + dx**2 / diffusivity_m2_s
            * (generation[1:-1] - recombination[1:-1])
        )
        # The outward surface flux equals S*Delta at x=0.
        residual[0] = (
            excess[1]
            - excess[0]
            - surface_velocity_m_s * dx / diffusivity_m2_s * excess[0]
        )
        residual[-1] = excess[-1] - right_excess_m3
        return residual / concentration_scale

    result = root(
        normalized_residual,
        initial,
        method="hybr",
        options={"xtol": 1.0e-10, "maxfev": 5000},
    )
    residual_max = np.max(np.abs(normalized_residual(result.x)))
    negativity_tolerance = 1.0e-10 * concentration_scale
    if (
        not result.success
        or residual_max > 2.0e-6
        or np.min(result.x) < -negativity_tolerance
    ):
        raise RuntimeError("steady SRH finite-difference solve did not converge")
    # Remove roundoff-level negative values at a zero-excess contact.
    return np.maximum(result.x, 0.0)


def dynamic_center_kinetics(
    time_s: np.ndarray,
    equilibrium_electron_m3: float,
    equilibrium_hole_m3: float,
    center_density_m3: float,
    electron_capture_m3_s: float,
    hole_capture_m3_s: float,
    electron_auxiliary_m3: float,
    hole_auxiliary_m3: float,
    generation_m3_s: float,
    pulse_end_s: float,
) -> dict[str, np.ndarray]:
    """Integrate carrier densities and center occupation through a generation pulse.

    Parameters
    ----------
    time_s : numpy.ndarray
        Strictly increasing output times in seconds.
    equilibrium_electron_m3, equilibrium_hole_m3 : float
        Positive initial equilibrium carrier densities in ``m^-3``.
    center_density_m3 : float
        Positive center density ``N_t`` in ``m^-3``.
    electron_capture_m3_s, hole_capture_m3_s : float
        Positive capture coefficients ``alpha_n`` and ``alpha_p`` in
        ``m^3 s^-1``.
    electron_auxiliary_m3, hole_auxiliary_m3 : float
        Positive detailed-balance densities ``n_B`` and ``p_B`` in ``m^-3``.
    generation_m3_s : float
        Nonnegative pair-generation rate during the pulse in ``m^-3 s^-1``.
    pulse_end_s : float
        Time at which generation is switched off, in seconds.

    Returns
    -------
    dict of str to numpy.ndarray
        Time histories named ``electron_m3``, ``hole_m3``, ``occupation``,
        ``electron_rate_m3_s``, ``hole_rate_m3_s``, and
        ``quasisteady_rate_m3_s``.

    Raises
    ------
    ValueError
        If times do not increase or a physical parameter is invalid.
    RuntimeError
        If the stiff BDF integration fails.

    Notes
    -----
    The integrated equations are ``dn/dt=G-U_n``, ``dp/dt=G-U_p``, and
    ``N_t*dF_t/dt=U_n-U_p``. Their charge invariant is
    ``p-n-N_t*F_t`` for the acceptor-like center convention. The algebraic SRH
    rate is evaluated alongside the dynamic rates to expose center-occupation
    lag.
    """
    times = np.asarray(time_s, dtype=float)
    if times.ndim != 1 or times.size < 2 or np.any(np.diff(times) <= 0.0):
        raise ValueError("time_s must be a strictly increasing one-dimensional grid")
    positive_parameters = (
        equilibrium_electron_m3,
        equilibrium_hole_m3,
        center_density_m3,
        electron_capture_m3_s,
        hole_capture_m3_s,
        electron_auxiliary_m3,
        hole_auxiliary_m3,
    )
    if any(value <= 0.0 for value in positive_parameters):
        raise ValueError("densities and capture coefficients must be positive")
    if generation_m3_s < 0.0 or not times[0] <= pulse_end_s <= times[-1]:
        raise ValueError("generation must be nonnegative and pulse end inside time grid")

    denominator = (
        electron_capture_m3_s
        * (equilibrium_electron_m3 + electron_auxiliary_m3)
        + hole_capture_m3_s
        * (equilibrium_hole_m3 + hole_auxiliary_m3)
    )
    initial_occupation = (
        electron_capture_m3_s * equilibrium_electron_m3
        + hole_capture_m3_s * hole_auxiliary_m3
    ) / denominator

    def rates(electron: float, hole: float, occupation: float) -> tuple[float, float]:
        electron_rate = center_density_m3 * electron_capture_m3_s * (
            electron * (1.0 - occupation) - electron_auxiliary_m3 * occupation
        )
        hole_rate = center_density_m3 * hole_capture_m3_s * (
            hole * occupation - hole_auxiliary_m3 * (1.0 - occupation)
        )
        return electron_rate, hole_rate

    def right_hand_side(time: float, state: np.ndarray) -> np.ndarray:
        electron, hole, occupation = state
        electron_rate, hole_rate = rates(electron, hole, occupation)
        generation = generation_m3_s if time < pulse_end_s else 0.0
        return np.array(
            (
                generation - electron_rate,
                generation - hole_rate,
                (electron_rate - hole_rate) / center_density_m3,
            )
        )

    integration = solve_ivp(
        right_hand_side,
        (times[0], times[-1]),
        (equilibrium_electron_m3, equilibrium_hole_m3, initial_occupation),
        method="BDF",
        t_eval=times,
        rtol=2.0e-8,
        atol=(1.0e4, 1.0e4, 1.0e-10),
        max_step=np.min(np.diff(times)),
    )
    if not integration.success:
        raise RuntimeError(f"dynamic center integration failed: {integration.message}")

    electron, hole, occupation = integration.y
    electron_rate, hole_rate = rates(electron, hole, occupation)
    intrinsic_squared = electron_auxiliary_m3 * hole_auxiliary_m3
    intrinsic = np.sqrt(intrinsic_squared)
    electron_lifetime = 1.0 / (center_density_m3 * electron_capture_m3_s)
    hole_lifetime = 1.0 / (center_density_m3 * hole_capture_m3_s)
    quasisteady = srh_recombination_rate(
        electron,
        hole,
        intrinsic,
        electron_lifetime,
        hole_lifetime,
        electron_auxiliary_m3,
        hole_auxiliary_m3,
    )
    return {
        "electron_m3": electron,
        "hole_m3": hole,
        "occupation": occupation,
        "electron_rate_m3_s": electron_rate,
        "hole_rate_m3_s": hole_rate,
        "quasisteady_rate_m3_s": quasisteady,
    }


def mobility_hall_figure() -> None:
    """Write the mobility, resistivity, and Hall-response figure.

    Returns
    -------
    None
        Writes ``transport-mobility-hall.pdf`` beneath :data:`FIGURE_DIR`.

    Notes
    -----
    The three panels compare Matthiessen mobility limits, the competing density
    and mobility contributions to resistivity, and the mobility-squared
    weighting of a weak-field two-carrier Hall coefficient.
    """
    temperature = np.linspace(100.0, 700.0, 500)
    impurity_levels = (1.0e15, 1.0e17, 1.0e19)  # cm^-3

    figure, axes = plt.subplots(1, 3, figsize=(11.0, 3.25), layout="constrained")

    # Panel 1: expose the opposite temperature trends of the two scattering laws.
    for impurity in impurity_levels:
        lattice = 1450.0 * (temperature / 300.0) ** -1.5
        ionized = 2200.0 * (temperature / 300.0) ** 1.5 * (1.0e16 / impurity)
        total = matthiessen_mobility(lattice, ionized)
        axes[0].loglog(temperature, total, label=fr"$N_I={impurity:.0e}\,\mathrm{{cm^{{-3}}}}$")
    axes[0].set(xlabel="temperature (K)", ylabel=r"mobility (cm$^2$ V$^{-1}$ s$^{-1}$)", title="Competing scattering limits")
    axes[0].legend(fontsize=7)
    axes[0].grid(alpha=0.2, which="both")

    # Panel 2: increasing density lowers mobility as it raises conductivity.
    donor = np.logspace(13.0, 20.0, 500)
    mobility_300 = matthiessen_mobility(
        np.full_like(donor, 1450.0),
        2200.0 * 1.0e16 / donor,
    )
    resistivity_ohm_cm = 1.0 / (Q * donor * mobility_300)
    axes[1].loglog(donor, resistivity_ohm_cm)
    axes[1].set(xlabel=r"electron density (cm$^{-3}$)", ylabel=r"resistivity ($\Omega$ cm)", title="Carrier density and mobility both matter")
    axes[1].grid(alpha=0.2, which="both")

    # Panel 3: normalize by n so only density and mobility ratios remain.
    density_ratio = np.logspace(-2.0, 2.0, 500)  # p/n
    electron_density = np.full_like(density_ratio, 1.0e22)
    hole_density = density_ratio * electron_density
    for mobility_ratio in (0.5, 1.0, 2.0):
        coefficient = two_carrier_hall_coefficient(
            electron_density,
            hole_density,
            electron_mobility=1.0,
            hole_mobility=mobility_ratio,
        )
        normalized = Q * electron_density * coefficient
        axes[2].semilogx(
            density_ratio,
            normalized,
            label=fr"$\mu_p/\mu_n={mobility_ratio:g}$",
        )
    axes[2].axhline(0.0, color="black", linewidth=0.7)
    axes[2].axvline(1.0, color="0.5", linewidth=0.7, linestyle=":")
    axes[2].set(
        xlabel=r"density ratio $p/n$",
        ylabel=r"normalized Hall response $qnR_H$",
        title="Two-carrier Hall sign",
    )
    axes[2].legend(fontsize=7)
    axes[2].grid(alpha=0.2)
    save(figure, "transport-mobility-hall.pdf")


def nonequilibrium_carriers_figure() -> None:
    """Write the nonequilibrium-carrier teaching figure.

    Returns
    -------
    None
        Writes ``transport-nonequilibrium-carriers.pdf`` beneath
        :data:`FIGURE_DIR`.

    Notes
    -----
    The panels show exponential lifetime decay, quasi-Fermi-level splitting,
    the energy dependence of Shockley--Read--Hall recombination, and
    Haynes--Shockley packet drift and diffusion.
    """
    temperature = 300.0
    intrinsic = 9.65e15
    majority = 1.0e22
    minority = intrinsic**2 / majority
    lifetime = 2.0e-6

    figure, axes = plt.subplots(2, 2, figsize=(9.2, 6.2), layout="constrained")
    # A pulsed low-injection experiment gives direct exponential lifetime decay.
    time = np.linspace(0.0, 8.0e-6, 500)
    initial = 2.0e20
    axes[0, 0].semilogy(time * 1.0e6, excess_decay(time, initial, lifetime) / initial)
    axes[0, 0].axvline(lifetime * 1.0e6, linestyle="--", color="0.3")
    axes[0, 0].set(xlabel="time after pulse (µs)", ylabel=r"$\Delta p/\Delta p_0$", title="Minority-carrier lifetime")
    axes[0, 0].grid(alpha=0.2)

    # Equal electron-hole injection separates the two quasi-Fermi levels.
    excess = np.logspace(12.0, 23.0, 500)
    electron = majority + excess
    hole = minority + excess
    splitting = quasi_fermi_splitting(electron, hole, intrinsic, temperature)
    axes[0, 1].semilogx(excess, splitting)
    axes[0, 1].set(xlabel=r"excess pairs $\Delta n=\Delta p$ (m$^{-3}$)", ylabel=r"$E_{Fn}-E_{Fp}$ (eV)", title="Quasi-Fermi-level splitting")
    axes[0, 1].grid(alpha=0.2)

    # The reduced SRH denominator is smallest for a trap near midgap.
    thermal_voltage = K_B * temperature / Q
    trap_offset = np.linspace(-0.45, 0.45, 600)
    n = majority + 1.0e20
    p = minority + 1.0e20
    srh_shape = (n * p - intrinsic**2) / (
        n + p + 2.0 * intrinsic * np.cosh(trap_offset / thermal_voltage)
    )
    srh_shape /= np.max(srh_shape)
    axes[1, 0].plot(trap_offset, srh_shape)
    axes[1, 0].set(
        xlabel=r"center energy $E_t-E_i$ (eV)",
        ylabel="normalized SRH rate",
        title="Equal-capture SRH limit",
    )
    axes[1, 0].grid(alpha=0.2)

    # The packet center drifts as v_d*t while its variance grows as 2*D*t.
    position = np.linspace(-0.02, 0.10, 700)
    diffusivity = 1.0e-3
    drift_velocity = 180.0
    pulse_times = (20.0e-6, 60.0e-6, 120.0e-6)
    for pulse_time in pulse_times:
        packet = np.exp(
            -(position - drift_velocity * pulse_time) ** 2
            / (4.0 * diffusivity * pulse_time)
            - pulse_time / (180.0e-6)
        ) / np.sqrt(4.0 * np.pi * diffusivity * pulse_time)
        packet /= np.max(packet)
        axes[1, 1].plot(position * 1.0e3, packet, label=fr"$t={pulse_time*1e6:.0f}\,\mu$s")
    axes[1, 1].set(xlabel="position (mm)", ylabel="normalized excess density", title="Haynes--Shockley drift and spreading")
    axes[1, 1].legend()
    axes[1, 1].grid(alpha=0.2)
    save(figure, "transport-nonequilibrium-carriers.pdf")


def continuity_numerics_figure() -> None:
    """Write finite-difference and dynamic-center validation results.

    Returns
    -------
    None
        Writes ``transport-continuity-numerics.pdf`` beneath
        :data:`FIGURE_DIR`.

    Notes
    -----
    The panels compare a finite-difference packet with its analytic
    Haynes--Shockley solution, solve a nonlinear steady SRH profile with surface
    recombination, and compare dynamic center rates with the algebraic
    steady-occupancy approximation.
    """
    figure, axes = plt.subplots(1, 3, figsize=(11.0, 3.35), layout="constrained")

    # Panel 1: use an infinite-domain analytic packet as a numerical oracle.
    position = np.linspace(-0.012, 0.012, 481, endpoint=False)
    times = np.linspace(0.0, 40.0e-6, 81)
    diffusivity = 8.0e-4
    drift_velocity = 150.0
    lifetime = 80.0e-6
    center = -4.0e-3
    initial_width = 0.8e-3
    initial = np.exp(-0.5 * ((position - center) / initial_width) ** 2)
    numerical = transient_continuity_fd(
        position,
        times,
        initial,
        diffusivity,
        drift_velocity,
        lifetime,
    )
    elapsed = times[-1]
    variance = initial_width**2 + 2.0 * diffusivity * elapsed
    analytic = (
        initial_width
        / np.sqrt(variance)
        * np.exp(
            -0.5
            * (position - center - drift_velocity * elapsed) ** 2
            / variance
            - elapsed / lifetime
        )
    )
    relative_l2 = np.linalg.norm(numerical[-1] - analytic) / np.linalg.norm(analytic)
    axes[0].plot(position * 1.0e3, analytic, label="analytic")
    axes[0].plot(position * 1.0e3, numerical[-1], "--", label="finite difference")
    axes[0].set(
        xlabel="position (mm)",
        ylabel="normalized excess density",
        title=fr"Packet validation: $L_2={relative_l2:.2e}$",
    )
    axes[0].legend(fontsize=7)
    axes[0].grid(alpha=0.2)

    # Panel 2: Beer--Lambert generation competes with bulk and surface loss.
    slab_position = np.linspace(0.0, 200.0e-6, 181)
    generation = 5.0e25 * np.exp(-slab_position / (40.0e-6))
    intrinsic = 1.0e16
    equilibrium_electron = 1.0e21
    equilibrium_hole = intrinsic**2 / equilibrium_electron
    for surface_velocity in (0.0, 100.0, 1000.0):
        profile = steady_srh_profile_fd(
            slab_position,
            generation,
            diffusivity_m2_s=1.0e-3,
            equilibrium_electron_m3=equilibrium_electron,
            equilibrium_hole_m3=equilibrium_hole,
            intrinsic_m3=intrinsic,
            electron_lifetime_s=2.0e-6,
            hole_lifetime_s=2.0e-6,
            electron_auxiliary_m3=intrinsic,
            hole_auxiliary_m3=intrinsic,
            surface_velocity_m_s=surface_velocity,
        )
        axes[1].semilogy(
            slab_position * 1.0e6,
            np.maximum(profile, 1.0),
            label=fr"$S={surface_velocity:g}$ m/s",
        )
    axes[1].set(
        xlabel="depth (µm)",
        ylabel=r"excess density (m$^{-3}$)",
        title="Nonlinear SRH boundary problem",
    )
    axes[1].legend(fontsize=7)
    axes[1].grid(alpha=0.2, which="both")

    # Panel 3: a short pulse exposes finite center-occupation response time.
    kinetic_times = np.linspace(0.0, 8.0e-6, 401)
    kinetic = dynamic_center_kinetics(
        kinetic_times,
        equilibrium_electron_m3=1.0e21,
        equilibrium_hole_m3=1.0e11,
        center_density_m3=2.0e20,
        electron_capture_m3_s=5.0e-15,
        hole_capture_m3_s=2.0e-14,
        electron_auxiliary_m3=1.0e16,
        hole_auxiliary_m3=1.0e16,
        generation_m3_s=2.0e26,
        pulse_end_s=0.6e-6,
    )
    rate_scale = 1.0e26
    axes[2].plot(
        kinetic_times * 1.0e6,
        kinetic["electron_rate_m3_s"] / rate_scale,
        label=r"dynamic $U_n$",
    )
    axes[2].plot(
        kinetic_times * 1.0e6,
        kinetic["hole_rate_m3_s"] / rate_scale,
        label=r"dynamic $U_p$",
    )
    axes[2].plot(
        kinetic_times * 1.0e6,
        kinetic["quasisteady_rate_m3_s"] / rate_scale,
        "--",
        label="algebraic SRH",
    )
    axes[2].axvline(0.6, color="0.4", linestyle=":", linewidth=0.8)
    axes[2].set(
        xlabel="time (µs)",
        ylabel=r"rate ($10^{26}$ m$^{-3}$ s$^{-1}$)",
        title="Dynamic center occupancy",
    )
    axes[2].legend(fontsize=7)
    axes[2].grid(alpha=0.2)
    save(figure, "transport-continuity-numerics.pdf")


def barrier_space_charge_numerics_figure() -> None:
    """Write finite-difference barrier, Poisson, and space-charge validations.

    Returns
    -------
    None
        Writes ``transport-barrier-sclc-numerics.pdf`` beneath
        :data:`FIGURE_DIR`.

    Notes
    -----
    The panels compare finite-difference quantum transmission with the exact
    rectangular-barrier result, validate numerical space-charge field profiles,
    and recover the Mott--Gurney and constant-velocity current laws.
    """
    figure, axes = plt.subplots(1, 3, figsize=(11.0, 3.35), layout="constrained")

    # Panel 1: convergence of the open-boundary finite-difference Hamiltonian.
    barrier_ev = 0.8
    energy_ev = np.linspace(0.03, 0.79, 180)
    exact = rectangular_barrier_transmission(
        energy_ev * Q,
        barrier_ev * Q,
        1.0e-9,
        0.26 * M_0,
    )
    axes[0].semilogy(energy_ev, exact, color="black", label="exact")
    for points in (20, 80):
        numerical = fd_rectangular_barrier_transmission(
            energy_ev * Q,
            barrier_ev * Q,
            1.0e-9,
            0.26 * M_0,
            points,
        )
        axes[0].semilogy(
            energy_ev,
            numerical,
            "--",
            label=fr"FD: $N={points}$",
        )
    axes[0].set(
        xlabel="incident energy (eV)",
        ylabel="transmission",
        title="Open-boundary Schrödinger grid",
        ylim=(1.0e-8, 1.2),
    )
    axes[0].legend(fontsize=7)
    axes[0].grid(alpha=0.2, which="both")

    # Panel 2: both transport limits satisfy the same voltage integral.
    voltage = 2.0
    length = 2.0e-6
    permittivity = 11.7 * EPSILON_0
    position_mu, field_mu, _ = space_charge_field_fd(
        voltage,
        length,
        permittivity,
        81,
        mobility_m2_v_s=0.12,
    )
    position_sat, field_sat, _ = space_charge_field_fd(
        voltage,
        length,
        permittivity,
        81,
        saturation_velocity_m_s=1.0e5,
    )
    reduced_position = position_mu / length
    axes[1].plot(
        reduced_position,
        field_mu / (voltage / length),
        label=r"constant $\mu$",
    )
    axes[1].plot(
        position_sat / length,
        field_sat / (voltage / length),
        label=r"constant $v_s$",
    )
    axes[1].plot(
        reduced_position,
        1.5 * np.sqrt(reduced_position),
        "k:",
        linewidth=1.0,
        label="analytic limits",
    )
    axes[1].plot(
        reduced_position,
        2.0 * reduced_position,
        "k:",
        linewidth=1.0,
    )
    axes[1].set(
        xlabel=r"reduced position $x/L$",
        ylabel=r"reduced field $\mathcal{E}/(V/L)$",
        title="Self-consistent space-charge field",
    )
    axes[1].legend(fontsize=7)
    axes[1].grid(alpha=0.2)

    # Panel 3: numerical currents recover the quadratic and linear voltage laws.
    voltage_sweep = np.logspace(-1.0, 1.0, 15)
    current_mu = []
    current_sat = []
    for applied_voltage in voltage_sweep:
        _, _, value_mu = space_charge_field_fd(
            applied_voltage,
            length,
            permittivity,
            41,
            mobility_m2_v_s=0.12,
        )
        _, _, value_sat = space_charge_field_fd(
            applied_voltage,
            length,
            permittivity,
            41,
            saturation_velocity_m_s=1.0e5,
        )
        current_mu.append(value_mu)
        current_sat.append(value_sat)
    axes[2].loglog(
        voltage_sweep,
        mott_gurney_current(voltage_sweep, permittivity, 0.12, length),
        color="C0",
        label=r"analytic $V^2$",
    )
    axes[2].loglog(
        voltage_sweep,
        current_mu,
        "o",
        ms=3,
        color="C0",
        label=r"FD constant $\mu$",
    )
    axes[2].loglog(
        voltage_sweep,
        saturated_space_charge_current(voltage_sweep, permittivity, 1.0e5, length),
        color="C1",
        label=r"analytic $V$",
    )
    axes[2].loglog(voltage_sweep, current_sat, "s", ms=3, color="C1", label=r"FD constant $v_s$")
    axes[2].set(
        xlabel="voltage (V)",
        ylabel=r"current density (A m$^{-2}$)",
        title="Finite-difference current laws",
    )
    axes[2].legend(fontsize=6.5)
    axes[2].grid(alpha=0.2, which="both")
    save(figure, "transport-barrier-sclc-numerics.pdf")


def barrier_space_charge_figure() -> None:
    """Write the barrier and space-charge transport figure.

    Returns
    -------
    None
        Writes ``transport-barriers-space-charge.pdf`` beneath
        :data:`FIGURE_DIR`.

    Notes
    -----
    The panels contrast thermal activation, coherent rectangular-barrier
    tunneling, and constant-mobility versus velocity-saturated
    space-charge-limited current.
    """
    temperature = 300.0
    thermal_voltage = K_B * temperature / Q
    figure, axes = plt.subplots(1, 3, figsize=(11.0, 3.25), layout="constrained")

    # Thermionic emission samples the exponentially small thermal population.
    barrier_ev = np.linspace(0.05, 1.2, 500)
    normalized_population = np.exp(-barrier_ev / thermal_voltage)
    axes[0].semilogy(barrier_ev, normalized_population)
    axes[0].set(xlabel="emission barrier (eV)", ylabel="Boltzmann tail fraction", title="Thermionic activation")
    axes[0].grid(alpha=0.2, which="both")

    # Coherent tunneling is exponentially sensitive to barrier width.
    barrier_ev_value = 0.8
    energy_ev = np.linspace(0.01, 0.79, 600)
    for width_nm in (0.6, 1.0, 1.6):
        transmission = rectangular_barrier_transmission(
            energy_ev * Q,
            barrier_ev_value * Q,
            width_nm * 1.0e-9,
            0.26 * M_0,
        )
        axes[1].semilogy(energy_ev, transmission, label=fr"$d={width_nm}$ nm")
    axes[1].set(xlabel="electron energy (eV)", ylabel="transmission", title="Rectangular-barrier tunneling", ylim=(1.0e-12, 1.2))
    axes[1].legend()
    axes[1].grid(alpha=0.2, which="both")

    # Space charge changes the current law from ohmic to injection limited.
    voltage = np.logspace(-2.0, 2.0, 500)
    permittivity = 11.7 * EPSILON_0
    mobility = 0.12
    length = 2.0e-6
    saturation_velocity = 1.0e5
    current_mg = mott_gurney_current(voltage, permittivity, mobility, length)
    current_sat = saturated_space_charge_current(voltage, permittivity, saturation_velocity, length)
    axes[2].loglog(voltage, current_mg, label=r"constant $\mu$: $J\propto V^2$")
    axes[2].loglog(voltage, current_sat, label=r"velocity saturation: $J\propto V$")
    axes[2].set(xlabel="voltage (V)", ylabel=r"current density (A m$^{-2}$)", title="Space-charge-limited current")
    axes[2].legend(fontsize=7)
    axes[2].grid(alpha=0.2, which="both")
    save(figure, "transport-barriers-space-charge.pdf")


def high_field_figure() -> None:
    """Write the high-field carrier-transport figure.

    Returns
    -------
    None
        Writes ``transport-high-field.pdf`` beneath :data:`FIGURE_DIR`.

    Notes
    -----
    The panels show ordinary velocity saturation, a reduced transferred-electron
    model with negative differential mobility, and empirical electron and hole
    impact-ionization coefficients.
    """
    figure, axes = plt.subplots(1, 3, figsize=(11.0, 3.25), layout="constrained")
    field = np.logspace(3.0, 8.0, 700)

    # Silicon electron and hole speeds approach different high-field limits.
    electron_velocity = saturation_velocity_model(field, 0.14, 1.0e5, 2.0)
    hole_velocity = saturation_velocity_model(field, 0.045, 0.8e5, 1.0)
    axes[0].loglog(field, electron_velocity, label="electrons")
    axes[0].loglog(field, hole_velocity, label="holes")
    axes[0].set(xlabel=r"field (V m$^{-1}$)", ylabel=r"drift speed (m s$^{-1}$)", title="Velocity saturation in silicon")
    axes[0].legend()
    axes[0].grid(alpha=0.2, which="both")

    # Transfer into a low-mobility satellite valley creates a falling branch.
    gaas_field = np.logspace(3.0, 7.2, 700)
    transfer_fraction = 1.0 / (1.0 + np.exp(-(np.log10(gaas_field) - 5.65) / 0.12))
    low_valley = saturation_velocity_model(gaas_field, 0.85, 2.1e5, 2.0)
    satellite = saturation_velocity_model(gaas_field, 0.08, 0.85e5, 1.0)
    gaas_velocity = (1.0 - transfer_fraction) * low_valley + transfer_fraction * satellite
    axes[1].loglog(gaas_field, gaas_velocity)
    axes[1].axvspan(3.5e5, 1.4e6, alpha=0.15, label="negative differential mobility")
    axes[1].set(xlabel=r"field (V m$^{-1}$)", ylabel=r"drift speed (m s$^{-1}$)", title="Transferred-electron effect")
    axes[1].legend(fontsize=7)
    axes[1].grid(alpha=0.2, which="both")

    # Empirical ionization coefficients activate rapidly above threshold fields.
    impact_field = np.linspace(1.0e7, 8.0e7, 600)
    alpha_n = ionization_rate(impact_field, 7.0e7, 1.23e8)
    alpha_p = ionization_rate(impact_field, 1.6e8, 2.04e8)
    axes[2].semilogy(impact_field / 1.0e7, alpha_n / 100.0, label=r"$\alpha_n$")
    axes[2].semilogy(impact_field / 1.0e7, alpha_p / 100.0, label=r"$\alpha_p$")
    axes[2].set(xlabel=r"field ($10^5$ V cm$^{-1}$)", ylabel=r"ionization rate (cm$^{-1}$)", title="Impact ionization")
    axes[2].legend()
    axes[2].grid(alpha=0.2, which="both")
    save(figure, "transport-high-field.pdf")
