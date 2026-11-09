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
    """Write one deterministic PDF beneath :data:`FIGURE_DIR` and close it."""
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
    """Combine two positive mobility limits with Matthiessen's rule.

    Parameters are expressed in the same mobility unit.  The returned mobility
    has that unit and is smaller than either individual limit.
    """
    lattice = np.asarray(lattice_mobility, dtype=float)
    impurity = np.asarray(impurity_mobility, dtype=float)
    if np.any(lattice <= 0.0) or np.any(impurity <= 0.0):
        raise ValueError("mobility limits must be positive")
    return 1.0 / (1.0 / lattice + 1.0 / impurity)


def hall_coefficient(density_m3: np.ndarray | float, charge_sign: int) -> np.ndarray:
    """Return the one-carrier Hall coefficient in cubic meters per coulomb."""
    density = np.asarray(density_m3, dtype=float)
    if np.any(density <= 0.0):
        raise ValueError("carrier density must be positive")
    if charge_sign not in (-1, 1):
        raise ValueError("charge_sign must be -1 for electrons or +1 for holes")
    return charge_sign / (Q * density)


def excess_decay(
    time_s: np.ndarray | float,
    initial_excess_m3: float,
    lifetime_s: float,
) -> np.ndarray:
    """Return low-injection excess-carrier decay after generation stops."""
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
    """Return ``E_Fn-E_Fp`` in electron-volts for nondegenerate carriers."""
    electron = np.asarray(electron_m3, dtype=float)
    hole = np.asarray(hole_m3, dtype=float)
    if np.any(electron <= 0.0) or np.any(hole <= 0.0):
        raise ValueError("carrier densities must be positive")
    if intrinsic_m3 <= 0.0 or temperature_k <= 0.0:
        raise ValueError("intrinsic density and temperature must be positive")
    thermal_voltage = K_B * temperature_k / Q
    return thermal_voltage * np.log(electron * hole / intrinsic_m3**2)


def rectangular_barrier_transmission(
    energy_j: np.ndarray | float,
    barrier_j: float,
    width_m: float,
    effective_mass_kg: float,
) -> np.ndarray:
    """Return the exact one-dimensional transmission for ``0 < E < V0``."""
    energy = np.asarray(energy_j, dtype=float)
    if np.any(energy <= 0.0) or np.any(energy >= barrier_j):
        raise ValueError("energy must satisfy 0 < energy < barrier")
    if barrier_j <= 0.0 or width_m <= 0.0 or effective_mass_kg <= 0.0:
        raise ValueError("barrier, width, and effective mass must be positive")
    kappa = np.sqrt(2.0 * effective_mass_kg * (barrier_j - energy)) / HBAR
    denominator = 1.0 + (
        barrier_j**2
        * np.sinh(np.minimum(kappa * width_m, 350.0)) ** 2
        / (4.0 * energy * (barrier_j - energy))
    )
    return 1.0 / denominator


def mott_gurney_current(
    voltage_v: np.ndarray | float,
    permittivity_f_m: float,
    mobility_m2_v_s: float,
    length_m: float,
) -> np.ndarray:
    """Return trap-free constant-mobility space-charge-limited current density."""
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
    """Return velocity-saturated space-charge-limited current density."""
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
    """Interpolate between linear drift and velocity saturation."""
    field = np.asarray(field_v_m, dtype=float)
    if np.any(field < 0.0):
        raise ValueError("field must be nonnegative")
    if mobility_m2_v_s <= 0.0 or saturation_velocity_m_s <= 0.0 or exponent <= 0.0:
        raise ValueError("mobility, saturation velocity, and exponent must be positive")
    reduced = mobility_m2_v_s * field / saturation_velocity_m_s
    return mobility_m2_v_s * field / (1.0 + reduced**exponent) ** (1.0 / exponent)


def ionization_rate(
    field_v_m: np.ndarray | float,
    prefactor_m1: float,
    characteristic_field_v_m: float,
) -> np.ndarray:
    """Return the empirical impact-ionization law ``alpha=A exp(-B/E)``."""
    field = np.asarray(field_v_m, dtype=float)
    if np.any(field <= 0.0):
        raise ValueError("field must be positive")
    if prefactor_m1 <= 0.0 or characteristic_field_v_m <= 0.0:
        raise ValueError("ionization parameters must be positive")
    return prefactor_m1 * np.exp(-characteristic_field_v_m / field)


def mobility_hall_figure() -> None:
    """Write mobility, resistivity, and Hall-response scaling curves."""
    temperature = np.linspace(100.0, 700.0, 500)
    impurity_levels = (1.0e15, 1.0e17, 1.0e19)  # cm^-3

    figure, axes = plt.subplots(1, 3, figsize=(11.0, 3.25), layout="constrained")
    for impurity in impurity_levels:
        lattice = 1450.0 * (temperature / 300.0) ** -1.5
        ionized = 2200.0 * (temperature / 300.0) ** 1.5 * (1.0e16 / impurity)
        total = matthiessen_mobility(lattice, ionized)
        axes[0].loglog(temperature, total, label=fr"$N_I={impurity:.0e}\,\mathrm{{cm^{{-3}}}}$")
    axes[0].set(xlabel="temperature (K)", ylabel=r"mobility (cm$^2$ V$^{-1}$ s$^{-1}$)", title="Competing scattering limits")
    axes[0].legend(fontsize=7)
    axes[0].grid(alpha=0.2, which="both")

    donor = np.logspace(13.0, 20.0, 500)
    mobility_300 = matthiessen_mobility(
        np.full_like(donor, 1450.0),
        2200.0 * 1.0e16 / donor,
    )
    resistivity_ohm_cm = 1.0 / (Q * donor * mobility_300)
    axes[1].loglog(donor, resistivity_ohm_cm)
    axes[1].set(xlabel=r"electron density (cm$^{-3}$)", ylabel=r"resistivity ($\Omega$ cm)", title="Carrier density and mobility both matter")
    axes[1].grid(alpha=0.2, which="both")

    density = np.logspace(20.0, 25.0, 400)
    current_density = 1.0e4
    magnetic_field = 0.5
    electron_field = hall_coefficient(density, -1) * current_density * magnetic_field
    hole_field = hall_coefficient(density, 1) * current_density * magnetic_field
    axes[2].semilogx(density, electron_field, label="electrons")
    axes[2].semilogx(density, hole_field, label="holes")
    axes[2].axhline(0.0, color="black", linewidth=0.7)
    axes[2].set(xlabel=r"carrier density (m$^{-3}$)", ylabel=r"Hall field $E_y$ (V m$^{-1}$)", title="Hall sign identifies carrier type")
    axes[2].legend()
    axes[2].grid(alpha=0.2)
    save(figure, "transport-mobility-hall.pdf")


def nonequilibrium_carriers_figure() -> None:
    """Write lifetime, quasi-Fermi, SRH, and Haynes--Shockley curves."""
    temperature = 300.0
    intrinsic = 9.65e15
    majority = 1.0e22
    minority = intrinsic**2 / majority
    lifetime = 2.0e-6

    figure, axes = plt.subplots(2, 2, figsize=(9.2, 6.2), layout="constrained")
    time = np.linspace(0.0, 8.0e-6, 500)
    initial = 2.0e20
    axes[0, 0].semilogy(time * 1.0e6, excess_decay(time, initial, lifetime) / initial)
    axes[0, 0].axvline(lifetime * 1.0e6, linestyle="--", color="0.3")
    axes[0, 0].set(xlabel="time after pulse (µs)", ylabel=r"$\Delta p/\Delta p_0$", title="Minority-carrier lifetime")
    axes[0, 0].grid(alpha=0.2)

    excess = np.logspace(12.0, 23.0, 500)
    electron = majority + excess
    hole = minority + excess
    splitting = quasi_fermi_splitting(electron, hole, intrinsic, temperature)
    axes[0, 1].semilogx(excess, splitting)
    axes[0, 1].set(xlabel=r"excess pairs $\Delta n=\Delta p$ (m$^{-3}$)", ylabel=r"$E_{Fn}-E_{Fp}$ (eV)", title="Quasi-Fermi-level splitting")
    axes[0, 1].grid(alpha=0.2)

    thermal_voltage = K_B * temperature / Q
    trap_offset = np.linspace(-0.45, 0.45, 600)
    n = majority + 1.0e20
    p = minority + 1.0e20
    srh_shape = (n * p - intrinsic**2) / (
        n + p + 2.0 * intrinsic * np.cosh(trap_offset / thermal_voltage)
    )
    srh_shape /= np.max(srh_shape)
    axes[1, 0].plot(trap_offset, srh_shape)
    axes[1, 0].set(xlabel=r"trap energy $E_t-E_i$ (eV)", ylabel="normalized SRH rate", title="Midgap centers dominate")
    axes[1, 0].grid(alpha=0.2)

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


def barrier_space_charge_figure() -> None:
    """Write thermionic, tunneling, and space-charge transport curves."""
    temperature = 300.0
    thermal_voltage = K_B * temperature / Q
    figure, axes = plt.subplots(1, 3, figsize=(11.0, 3.25), layout="constrained")

    barrier_ev = np.linspace(0.05, 1.2, 500)
    normalized_population = np.exp(-barrier_ev / thermal_voltage)
    axes[0].semilogy(barrier_ev, normalized_population)
    axes[0].set(xlabel="emission barrier (eV)", ylabel="Boltzmann tail fraction", title="Thermionic activation")
    axes[0].grid(alpha=0.2, which="both")

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
    """Write velocity saturation, transferred-electron, and ionization curves."""
    figure, axes = plt.subplots(1, 3, figsize=(11.0, 3.25), layout="constrained")
    field = np.logspace(3.0, 8.0, 700)

    electron_velocity = saturation_velocity_model(field, 0.14, 1.0e5, 2.0)
    hole_velocity = saturation_velocity_model(field, 0.045, 0.8e5, 1.0)
    axes[0].loglog(field, electron_velocity, label="electrons")
    axes[0].loglog(field, hole_velocity, label="holes")
    axes[0].set(xlabel=r"field (V m$^{-1}$)", ylabel=r"drift speed (m s$^{-1}$)", title="Velocity saturation in silicon")
    axes[0].legend()
    axes[0].grid(alpha=0.2, which="both")

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

    impact_field = np.linspace(1.0e7, 8.0e7, 600)
    alpha_n = ionization_rate(impact_field, 7.0e7, 1.23e8)
    alpha_p = ionization_rate(impact_field, 1.6e8, 2.04e8)
    axes[2].semilogy(impact_field / 1.0e7, alpha_n / 100.0, label=r"$\alpha_n$")
    axes[2].semilogy(impact_field / 1.0e7, alpha_p / 100.0, label=r"$\alpha_p$")
    axes[2].set(xlabel=r"field ($10^5$ V cm$^{-1}$)", ylabel=r"ionization rate (cm$^{-1}$)", title="Impact ionization")
    axes[2].legend()
    axes[2].grid(alpha=0.2, which="both")
    save(figure, "transport-high-field.pdf")
