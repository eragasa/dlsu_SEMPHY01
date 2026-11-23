"""Integration tools connecting equilibrium states to transport and junctions.

This module supports Chapter 21 of the SEMPHY01 notes.  It makes the data
contract between equilibrium carrier statistics and device prediction explicit,
including the generalized Einstein relation, chemical-potential built-in
voltage, ionized-dopant depletion charge, and uncertainty ranking.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import Boltzmann as K_B
from scipy.constants import elementary_charge as Q
from scipy.integrate import quad
from scipy.special import expit, gamma

from .junctions import FIGURE_DIR

mpl.rcParams.update(
    {
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "legend.fontsize": 8,
        "figure.dpi": 140,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
    }
)


def built_in_from_chemical_potentials(
    n_bulk_chemical_potential_eV: float,
    p_bulk_chemical_potential_eV: float,
) -> float:
    """Calculate homojunction built-in voltage from isolated bulk states.

    Parameters
    ----------
    n_bulk_chemical_potential_eV : float
        Chemical potential of the isolated n-type bulk in electron-volts.
    p_bulk_chemical_potential_eV : float
        Chemical potential of the isolated p-type bulk in electron-volts,
        using the same energy reference.

    Returns
    -------
    float
        Built-in voltage in volts.

    Notes
    -----
    An energy difference of one electron-volt divided by the positive
    elementary charge is numerically one volt.  A common rigid shift of both
    chemical potentials therefore leaves the result invariant.
    """
    return float(n_bulk_chemical_potential_eV - p_bulk_chemical_potential_eV)


def built_in_from_carrier_states(
    n_side_electron_m3: float,
    p_side_hole_m3: float,
    intrinsic_m3: float,
    temperature_K: float,
) -> float:
    """Calculate the nondegenerate built-in voltage from bulk carrier states.

    Parameters
    ----------
    n_side_electron_m3 : float
        Equilibrium majority-electron density on the n-side in ``m^-3``.
    p_side_hole_m3 : float
        Equilibrium majority-hole density on the p-side in ``m^-3``.
    intrinsic_m3 : float
        Intrinsic carrier density in ``m^-3``.
    temperature_K : float
        Positive absolute temperature in kelvin.

    Returns
    -------
    float
        Built-in voltage in volts.
    """
    if min(n_side_electron_m3, p_side_hole_m3, intrinsic_m3, temperature_K) <= 0.0:
        raise ValueError("carrier densities and temperature must be positive")
    thermal_voltage = K_B * temperature_K / Q
    return float(
        thermal_voltage
        * np.log(n_side_electron_m3 * p_side_hole_m3 / intrinsic_m3**2)
    )


def normalized_fermi_integral(order: float, reduced_potential: float) -> float:
    """Evaluate a normalized complete Fermi--Dirac integral.

    Parameters
    ----------
    order : float
        Integral order greater than ``-1``.
    reduced_potential : float
        Reduced chemical potential ``eta``.

    Returns
    -------
    float
        Normalized Fermi integral
        ``F_j(eta)=integral[x^j/(1+exp(x-eta))] dx / Gamma(j+1)``.
    """
    if order <= -1.0:
        raise ValueError("order must exceed -1")
    value, _ = quad(
        lambda energy: energy**order * expit(reduced_potential - energy),
        0.0,
        np.inf,
        epsabs=1.0e-11,
        epsrel=1.0e-10,
        limit=300,
    )
    return float(value / gamma(order + 1.0))


def generalized_einstein_factor(reduced_potential: np.ndarray | float) -> np.ndarray:
    """Calculate the degeneracy correction to the Einstein voltage.

    Parameters
    ----------
    reduced_potential : array_like
        Reduced electron chemical potential ``eta=(mu-E_C)/(k_B T)``.

    Returns
    -------
    numpy.ndarray
        Dimensionless factor ``F_(1/2)(eta)/F_(-1/2)(eta)`` such that
        ``D/mu=V_T*factor`` for a three-dimensional parabolic band.
    """
    eta = np.asarray(reduced_potential, dtype=float)
    flat = eta.ravel()
    factor = np.array(
        [
            normalized_fermi_integral(0.5, value)
            / normalized_fermi_integral(-0.5, value)
            for value in flat
        ]
    )
    return factor.reshape(eta.shape)


def depletion_width_from_ionized_dopants(
    built_in_voltage_V: float,
    applied_voltage_V: float,
    ionized_acceptor_m3: float,
    ionized_donor_m3: float,
    permittivity_F_m: float,
) -> tuple[float, float, float]:
    """Calculate depletion widths using transferred ionized dopant densities.

    Parameters
    ----------
    built_in_voltage_V : float
        Positive built-in potential in volts.
    applied_voltage_V : float
        Junction voltage in volts, positive for forward bias.
    ionized_acceptor_m3, ionized_donor_m3 : float
        Positive ionized acceptor and donor densities in ``m^-3``.
    permittivity_F_m : float
        Positive semiconductor permittivity in ``F m^-1``.

    Returns
    -------
    total_width_m : float
        Total depletion width in meters.
    p_width_m : float
        Depletion width on the p-side in meters.
    n_width_m : float
        Depletion width on the n-side in meters.
    """
    if min(
        built_in_voltage_V,
        ionized_acceptor_m3,
        ionized_donor_m3,
        permittivity_F_m,
    ) <= 0.0:
        raise ValueError("built-in voltage, ionized densities, and permittivity must be positive")
    drop = built_in_voltage_V - applied_voltage_V
    if drop <= 0.0:
        raise ValueError("applied forward voltage must remain below built-in voltage")
    total = np.sqrt(
        2.0
        * permittivity_F_m
        * drop
        / Q
        * (1.0 / ionized_acceptor_m3 + 1.0 / ionized_donor_m3)
    )
    p_width = ionized_donor_m3 / (ionized_acceptor_m3 + ionized_donor_m3) * total
    n_width = ionized_acceptor_m3 / (ionized_acceptor_m3 + ionized_donor_m3) * total
    return float(total), float(p_width), float(n_width)


def local_log_sensitivity(
    baseline: float,
    perturbed: float,
    fractional_parameter_change: float,
) -> float:
    """Estimate a local logarithmic output sensitivity by finite difference.

    Parameters
    ----------
    baseline : float
        Positive baseline output.
    perturbed : float
        Positive output after a parameter perturbation.
    fractional_parameter_change : float
        Nonzero fractional parameter perturbation, for example ``0.01``.

    Returns
    -------
    float
        Approximation to ``d ln(output) / d ln(parameter)``.
    """
    if baseline <= 0.0 or perturbed <= 0.0 or fractional_parameter_change == 0.0:
        raise ValueError("outputs must be positive and perturbation nonzero")
    return float(np.log(perturbed / baseline) / np.log1p(fractional_parameter_change))


def integrated_prediction_figure() -> None:
    """Write the Chapter 21 integration and uncertainty figure."""
    figure, axes = plt.subplots(2, 2, figsize=(10.5, 7.0), layout="constrained")

    eta = np.linspace(-8.0, 8.0, 100)
    factor = generalized_einstein_factor(eta)
    axes[0, 0].plot(eta, factor, color="#1769aa")
    axes[0, 0].axhline(1.0, color="black", linestyle="--", label="Boltzmann limit")
    axes[0, 0].set(
        xlabel=r"reduced potential $\eta=(\mu-E_C)/(k_BT)$",
        ylabel=r"$(D/\mu)/V_T$",
        title="Generalized Einstein relation",
    )
    axes[0, 0].legend()

    ionization_fraction = np.linspace(0.1, 1.0, 200)
    # A one-sided junction has W proportional to the inverse square root of the
    # lightly doped side's ionized fraction and C proportional to 1/W.
    axes[0, 1].plot(ionization_fraction, ionization_fraction ** -0.5, label=r"$W/W_{\rm full}$")
    axes[0, 1].plot(ionization_fraction, ionization_fraction ** 0.5, label=r"$C/C_{\rm full}$")
    axes[0, 1].set(
        xlabel="light-side ionized fraction",
        ylabel="ratio to complete-ionization result",
        title="Incomplete ionization changes depletion",
    )
    axes[0, 1].legend()

    temperature = 300.0
    thermal_voltage = K_B * temperature / Q
    intrinsic = 1.0e16
    n_side = 1.0e22
    p_side = 1.0e23
    mu_n = thermal_voltage * np.log(n_side / intrinsic)
    mu_p = -thermal_voltage * np.log(p_side / intrinsic)
    shifts = np.linspace(-2.0, 2.0, 21)
    chemical = [built_in_from_chemical_potentials(mu_n + shift, mu_p + shift) for shift in shifts]
    carrier = built_in_from_carrier_states(n_side, p_side, intrinsic, temperature)
    axes[1, 0].plot(shifts, chemical, "o", markersize=3, label="chemical-potential difference")
    axes[1, 0].axhline(carrier, color="#b23a48", linestyle="--", label="carrier-density expression")
    axes[1, 0].set(
        xlabel="common energy-reference shift (eV)",
        ylabel=r"$V_{bi}$ (V)",
        title="Physical voltage is energy-reference invariant",
    )
    axes[1, 0].legend()

    labels = [r"$E_g\pm10$ meV", r"$\mu\pm10\%$", r"$\tau\pm50\%$", r"$N\pm5\%$"]
    changes = [100.0 * 0.010 / thermal_voltage, 5.0, 25.0, 5.0]
    axes[1, 1].barh(labels, changes, color=["#b23a48", "#1769aa", "#237a57", "#d97706"])
    axes[1, 1].set(
        xlabel="approximate magnitude of relative change in ideal $J_s$ (%)",
        title="Declared uncertainties produce unequal effects",
    )

    for axis in axes.flat:
        axis.grid(alpha=0.2)
    figure.savefig(
        FIGURE_DIR / "device-prediction-integration.pdf",
        metadata={
            "Creator": "SEMPHY01 device-prediction integration notebook",
            "Producer": "Matplotlib",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    plt.close(figure)
