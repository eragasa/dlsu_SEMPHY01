"""Original schematic figure atlas for Sze--Lee Chapter 3 topics.

The figures are generated from equations or redrawn as original instructional
schematics.  They cover the physical subjects of Figs. 1--30 in Chapter 3 of
Sze and Lee without reproducing the textbook artwork.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Rectangle
from scipy.constants import Boltzmann as K_B
from scipy.constants import elementary_charge as Q
from scipy.constants import epsilon_0 as EPSILON_0
from scipy.optimize import brentq

from .junctions import (
    FIGURE_DIR,
    abrupt_depletion_profile,
    avalanche_breakdown_voltage,
    ideal_diffusion_saturation_current_density,
    ideal_diode_current_density,
    solve_poisson_boltzmann_junction,
)

mpl.rcParams.update(
    {
        "font.size": 8.5,
        "axes.labelsize": 8.5,
        "axes.titlesize": 9.5,
        "legend.fontsize": 7,
        "figure.dpi": 140,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
    }
)

BLUE = "#1769aa"
RED = "#b23a48"
GREEN = "#237a57"
ORANGE = "#d97706"
GRAY = "#5f6368"


def _save(figure: plt.Figure, name: str) -> None:
    """Write a deterministic PDF to the manuscript figure directory."""
    figure.savefig(
        FIGURE_DIR / name,
        metadata={
            "Creator": "SEMPHY01 original junction schematic generator",
            "Producer": "Matplotlib",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    plt.close(figure)


def _clean(axis: plt.Axes) -> None:
    """Remove upper and right spines and add a light grid."""
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(alpha=0.17)


def _band_diagram(axis: plt.Axes, barrier: float, title: str, *, offsets: bool = False) -> None:
    """Draw a compact p--n band diagram with a specified barrier height."""
    x = np.linspace(-1.0, 1.0, 300)
    bend = barrier * (1.0 - np.tanh(3.0 * x)) / 2.0
    if offsets:
        bend = bend + np.where(x >= 0.0, 0.18, 0.0)
    axis.plot(x, bend + 0.55, color=BLUE, label=r"$E_C$")
    axis.plot(x, bend - 0.55, color=RED, label=r"$E_V$")
    axis.axhline(0.0, color=GRAY, linestyle="--", label=r"$E_F$")
    axis.axvline(0.0, color="0.65", linestyle=":")
    axis.text(-0.86, -0.92, "p", color=RED, fontweight="bold")
    axis.text(0.80, -0.92, "n", color=BLUE, fontweight="bold")
    axis.set(xlim=(-1, 1), ylim=(-1.0, 1.45), xlabel="position", ylabel="energy", title=title)
    axis.set_yticks([])
    _clean(axis)


def junction_formation_atlas() -> None:
    """Generate one parameter-consistent junction-formation calculation.

    Every panel uses the same nondegenerate Si-like abrupt junction.  Isolated
    Fermi levels follow effective-density-of-states statistics, equilibrium
    bands and charge come from the nonlinear Poisson--Boltzmann solution, the
    comparison profile is the analytic depletion approximation, and the I--V
    curve uses the long-base minority-diffusion saturation current.
    """
    acceptor, donor, intrinsic, temperature = 5.0e22, 1.0e22, 1.0e16, 300.0
    permittivity = 11.7 * EPSILON_0
    thermal_voltage = K_B * temperature / Q
    conduction_dos, valence_dos = 2.8e25, 1.04e25
    # Enforce n_i^2=N_C*N_V*exp(-E_g/kT), so the isolated Fermi-level
    # separation and the Poisson--Boltzmann contact potential are identical.
    band_gap_ev = thermal_voltage * np.log(
        conduction_dos * valence_dos / intrinsic**2
    )
    depletion = abrupt_depletion_profile(
        acceptor, donor, intrinsic, temperature, permittivity
    )
    numerical = solve_poisson_boltzmann_junction(
        acceptor,
        donor,
        intrinsic,
        temperature,
        permittivity,
        points=601,
    )

    electron_mobility, hole_mobility = 0.135, 0.048
    electron_lifetime = hole_lifetime = 1.0e-6
    electron_diffusivity = electron_mobility * thermal_voltage
    hole_diffusivity = hole_mobility * thermal_voltage
    saturation_current = ideal_diffusion_saturation_current_density(
        acceptor,
        donor,
        intrinsic,
        electron_diffusivity,
        hole_diffusivity,
        electron_lifetime,
        hole_lifetime,
    )

    figure, axes = plt.subplots(2, 2, figsize=(9.0, 6.2), layout="constrained")
    voltage = np.linspace(-0.20, 0.72, 700)
    current = ideal_diode_current_density(voltage, saturation_current, temperature)
    axes[0, 0].plot(voltage, current, color=BLUE)
    axes[0, 0].axhline(-saturation_current, color=GRAY, linestyle="--", label=r"$-J_s$")
    axes[0, 0].axhline(0, color="black", linewidth=0.7)
    axes[0, 0].axvline(0, color="black", linewidth=0.7)
    axes[0, 0].set_yscale("symlog", linthresh=saturation_current / 10.0)
    axes[0, 0].set(
        xlabel="junction voltage (V)",
        ylabel=r"current density (A m$^{-2}$)",
        title="Long-base diffusion rectification",
    )
    axes[0, 0].legend()
    _clean(axes[0, 0])

    p_fermi_ev = thermal_voltage * np.log(valence_dos / acceptor)
    n_fermi_ev = band_gap_ev - thermal_voltage * np.log(conduction_dos / donor)
    left_x = np.linspace(-1.0, -0.03, 150)
    right_x = np.linspace(0.03, 1.0, 150)
    axes[0, 1].plot((-1.0, 1.0), (band_gap_ev, band_gap_ev), color=BLUE, label=r"$E_C$")
    axes[0, 1].plot((-1.0, 1.0), (0.0, 0.0), color=RED, label=r"$E_V$")
    axes[0, 1].plot(left_x, np.full_like(left_x, p_fermi_ev), color=GRAY, linewidth=1.8)
    axes[0, 1].plot(right_x, np.full_like(right_x, n_fermi_ev), color=GRAY, linewidth=1.8, label=r"$E_F$")
    axes[0, 1].axvline(0, color="0.6", linestyle=":")
    axes[0, 1].annotate(
        "",
        xy=(0.34, n_fermi_ev),
        xytext=(0.34, p_fermi_ev),
        arrowprops={"arrowstyle": "<->"},
    )
    axes[0, 1].text(
        0.39,
        0.5 * (p_fermi_ev + n_fermi_ev),
        fr"$qV_{{bi}}={float(depletion['built_in_v']):.3f}\,$eV",
        va="center",
    )
    axes[0, 1].text(-0.80, 0.05, "p", color=RED, fontweight="bold")
    axes[0, 1].text(0.72, 0.05, "n", color=BLUE, fontweight="bold")
    axes[0, 1].set(
        xlabel="isolated samples",
        ylabel="energy (eV)",
        title="Before contact: carrier-statistics levels",
        ylim=(-0.08, band_gap_ev + 0.10),
    )
    axes[0, 1].legend(ncol=3, loc="upper center")
    _clean(axes[0, 1])

    position_um = np.asarray(numerical["position_m"]) * 1.0e6
    potential = np.asarray(numerical["potential_v"])
    conduction_band = thermal_voltage * np.log(conduction_dos / intrinsic) - potential
    valence_band = conduction_band - band_gap_ev
    axes[1, 0].plot(position_um, conduction_band, color=BLUE, label=r"$E_C$")
    axes[1, 0].plot(position_um, valence_band, color=RED, label=r"$E_V$")
    axes[1, 0].axhline(0.0, color=GRAY, linestyle="--", label=r"$E_F$")
    axes[1, 0].axvline(0.0, color="0.65", linestyle=":")
    axes[1, 0].set(
        xlabel="position (µm)",
        ylabel="energy relative to $E_F$ (eV)",
        title="After contact: Poisson--Boltzmann bands",
    )
    axes[1, 0].legend(ncol=3)
    _clean(axes[1, 0])

    position_m = np.asarray(numerical["position_m"])
    depletion_charge = np.zeros_like(position_m)
    p_width = float(depletion["p_width_m"])
    n_width = float(depletion["n_width_m"])
    depletion_charge[(position_m >= -p_width) & (position_m < 0.0)] = -Q * acceptor
    depletion_charge[(position_m >= 0.0) & (position_m <= n_width)] = Q * donor
    axes[1, 1].plot(
        position_um,
        np.asarray(numerical["charge_c_m3"]) / 1.0e3,
        color=GREEN,
        label="Poisson--Boltzmann",
    )
    axes[1, 1].step(
        position_um,
        depletion_charge / 1.0e3,
        where="mid",
        color=ORANGE,
        linestyle="--",
        label="depletion approximation",
    )
    axes[1, 1].axhline(0, color="black", linewidth=0.7)
    axes[1, 1].set(
        xlabel="position (µm)",
        ylabel=r"charge density ($10^3$ C m$^{-3}$)",
        title="Calculated space-charge transition",
    )
    axes[1, 1].legend()
    _clean(axes[1, 1])
    figure.suptitle(
        r"One consistent 300 K junction: "
        r"$N_A=5\times10^{16}$ cm$^{-3}$, "
        r"$N_D=10^{16}$ cm$^{-3}$",
        fontsize=10,
    )
    _save(figure, "junction-atlas-formation.pdf")


def doping_electrostatics_atlas() -> None:
    """Generate built-in-potential, doping-profile, and bias-state panels."""
    figure, axes = plt.subplots(2, 2, figsize=(9.0, 6.2), layout="constrained")
    density_cm3 = np.logspace(14, 19, 300)
    thermal_voltage = K_B * 300.0 / Q
    for intrinsic_cm3, label, style in ((1.0e10, "Si", "-"), (2.0e6, "GaAs", "--")):
        component = thermal_voltage * np.log(density_cm3 / intrinsic_cm3)
        axes[0, 0].semilogx(density_cm3, component, style, label=label)
    axes[0, 0].set(xlabel=r"dopant density (cm$^{-3}$)", ylabel="bulk potential magnitude (V)", title="Neutral-region potential components")
    axes[0, 0].legend()
    _clean(axes[0, 0])

    x = np.linspace(-1.0, 1.0, 401)
    abrupt = np.where(x < 0.0, -1.0, 1.0)
    graded = x
    axes[0, 1].plot(x, abrupt, label="abrupt")
    axes[0, 1].plot(x, graded, label="linearly graded")
    axes[0, 1].axhline(0, color="black", linewidth=0.7)
    axes[0, 1].set(xlabel="position", ylabel=r"$N_D-N_A$ (normalized)", title="Canonical doping profiles")
    axes[0, 1].legend()
    _clean(axes[0, 1])

    for barrier, label, color in ((0.75, "equilibrium", GRAY), (0.38, "forward bias", RED), (1.05, "reverse bias", BLUE)):
        bend = barrier * (1.0 - np.tanh(3.5 * x)) / 2
        axes[1, 0].plot(x, bend, color=color, label=label)
    axes[1, 0].set(xlabel="position", ylabel=r"$E_C$ shift (eV)", title="Bias changes barrier and depletion width")
    axes[1, 0].legend()
    _clean(axes[1, 0])

    widths = (0.62, 0.42, 0.85)
    for width, label, color in zip(widths, ("equilibrium", "forward", "reverse"), (GRAY, RED, BLUE), strict=True):
        axes[1, 1].plot([-width, 0.0, width], [0.0, -1.0, 0.0], color=color, label=label)
    axes[1, 1].set(xlabel="position", ylabel="field (normalized)", title="Bias-dependent depletion field")
    axes[1, 1].legend()
    _clean(axes[1, 1])
    _save(figure, "junction-atlas-doping-bias.pdf")


def graded_junction_atlas() -> None:
    """Generate linearly graded charge, field, potential, and material trends."""
    figure, axes = plt.subplots(2, 2, figsize=(9.0, 6.2), layout="constrained")
    x = np.linspace(-1.0, 1.0, 501)
    charge = x
    field = 0.5 * (x**2 - 1.0)
    potential = 0.5 * x - x**3 / 6.0
    axes[0, 0].plot(x, charge, color=GREEN)
    axes[0, 0].axhline(0, color="black", linewidth=0.7)
    axes[0, 0].set(xlabel="position", ylabel="net doping", title="Linear impurity gradient")
    axes[0, 1].plot(x, field, color=BLUE)
    axes[0, 1].set(xlabel="position", ylabel="field", title="Parabolic electric field")
    axes[1, 0].plot(x, potential - potential.min(), color=ORANGE, label=r"$\psi$")
    axes[1, 0].plot(x, 0.7 - (potential - potential.min()), color=BLUE, label=r"$E_C$ (schematic)")
    axes[1, 0].set(xlabel="position", ylabel="normalized value", title="Cubic potential and band bending")
    axes[1, 0].legend()

    gradients_cm4 = np.logspace(18, 24, 100)
    for ni_cm3, eps_r, label in ((1e10, 11.7, "Si"), (2e6, 12.9, "GaAs")):
        values = []
        for gradient_cm4 in gradients_cm4:
            gradient_m4 = gradient_cm4 * 1e8
            ni_m3 = ni_cm3 * 1e6
            eps = eps_r * EPSILON_0
            def equation(voltage: float) -> float:
                width = (12.0 * eps * voltage / (Q * gradient_m4)) ** (1.0 / 3.0)
                edge_doping = gradient_m4 * width / 2.0
                return voltage - 2.0 * (K_B * 300.0 / Q) * np.log(edge_doping / ni_m3)
            values.append(brentq(equation, 0.05, 2.5))
        axes[1, 1].semilogx(gradients_cm4, values, label=label)
    axes[1, 1].set(xlabel=r"impurity gradient (cm$^{-4}$)", ylabel=r"$V_{bi}$ (V)", title="Graded-junction built-in potential")
    axes[1, 1].legend()
    for axis in axes.flat:
        _clean(axis)
    _save(figure, "junction-atlas-graded.pdf")


def capacitance_construction_atlas() -> None:
    """Generate incremental depletion-charge and profile-engineering panels."""
    figure, axes = plt.subplots(1, 3, figsize=(11.0, 3.4), layout="constrained")
    x = np.linspace(0.0, 1.2, 400)
    doping = 0.55 + 0.35 * x + 0.08 * np.sin(7 * x)
    w1, w2 = 0.72, 0.88
    axes[0].plot(x, doping, color=GREEN)
    axes[0].fill_between(x, 0, doping, where=(x >= w1) & (x <= w2), color=ORANGE, alpha=0.45, label=r"$\mathrm{d}Q/q$")
    axes[0].axvline(w1, color=GRAY, linestyle="--")
    axes[0].axvline(w2, color=GRAY, linestyle="--")
    axes[0].set(xlabel="light-side position", ylabel="dopant density", title="Moving depletion edge")
    axes[0].legend()

    field1 = -np.maximum(w1 - x, 0) * (0.55 + 0.15 * x)
    field2 = -np.maximum(w2 - x, 0) * (0.55 + 0.15 * x)
    axes[1].plot(x, field1, label=r"$V_R$")
    axes[1].plot(x, field2, label=r"$V_R+\mathrm{d}V$")
    axes[1].fill_between(x, field1, field2, color=ORANGE, alpha=0.35)
    axes[1].set(xlabel="position", ylabel="electric field", title="Incremental field area")
    axes[1].legend()

    centered = np.linspace(-1, 1, 400)
    profiles = {
        "hyperabrupt": np.sign(centered) * np.maximum(np.abs(centered), 0.08) ** -0.45,
        "abrupt": np.sign(centered),
        "linearly graded": centered,
    }
    for label, profile in profiles.items():
        axes[2].plot(centered, profile, label=label)
    axes[2].set_ylim(-3.3, 3.3)
    axes[2].set(xlabel="position", ylabel="net doping (normalized)", title="Varactor profile engineering")
    axes[2].legend()
    for axis in axes:
        _clean(axis)
    _save(figure, "junction-atlas-capacitance-construction.pdf")


def injection_transport_atlas() -> None:
    """Generate biased band, carrier-injection, and current-component panels."""
    figure, axes = plt.subplots(2, 2, figsize=(9.0, 6.2), layout="constrained")
    _band_diagram(axes[0, 0], 0.36, "Forward bias: reduced barrier")
    _band_diagram(axes[0, 1], 1.04, "Reverse bias: increased barrier")
    x = np.linspace(-3.0, 3.0, 601)
    dn = np.where(x < 0, np.exp(x / 0.75), 0.0)
    dp = np.where(x >= 0, np.exp(-x / 0.95), 0.0)
    axes[1, 0].plot(x, dn, color=BLUE, label=r"injected $\Delta n_p$")
    axes[1, 0].plot(x, dp, color=RED, label=r"injected $\Delta p_n$")
    axes[1, 0].fill_between(x, 0, dn + dp, color=GREEN, alpha=0.15, label="stored minority charge")
    axes[1, 0].axvspan(-0.18, 0.18, color="0.8", alpha=0.7, label="depletion region")
    axes[1, 0].set(xlabel="position", ylabel="excess density", title="Forward-bias minority injection")
    axes[1, 0].legend(ncol=2)

    electron = 0.5 * (1.0 + np.tanh(1.5 * x))
    hole = 1.0 - electron
    axes[1, 1].plot(x, electron, color=BLUE, label=r"$J_n$")
    axes[1, 1].plot(x, hole, color=RED, label=r"$J_p$")
    axes[1, 1].plot(x, electron + hole, color="black", linestyle="--", label=r"$J_n+J_p$")
    axes[1, 1].set(xlabel="position", ylabel="current component", title="Ideal steady current conversion")
    axes[1, 1].legend()
    for axis in axes.flat:
        _clean(axis)
    _save(figure, "junction-atlas-injection.pdf")


def diode_characteristics_atlas() -> None:
    """Generate ideal, nonideal, and temperature-dependent diode characteristics."""
    figure, axes = plt.subplots(2, 2, figsize=(9.2, 6.4), layout="constrained")
    thermal_voltage = K_B * 300.0 / Q
    voltage = np.linspace(-0.25, 0.82, 500)
    current = 1e-11 * np.expm1(np.clip(voltage / thermal_voltage, -100, 28))
    axes[0, 0].plot(voltage, current, color=BLUE)
    axes[0, 0].set(xlabel="voltage (V)", ylabel="current density", title="Ideal diode: Cartesian scale")
    axes[0, 0].set_ylim(-2e-11, 18)
    axes[0, 1].semilogy(voltage, np.maximum(np.abs(current), 1e-14), color=BLUE)
    axes[0, 1].set(xlabel="voltage (V)", ylabel="current magnitude", title="Ideal diode: semilog scale")

    forward = np.linspace(0.05, 1.25, 500)
    for material, js, j02, resistance, color in (("Si", 1e-11, 2e-8, 0.035, BLUE), ("GaAs", 1e-16, 2e-10, 0.06, RED)):
        ideal = js * np.exp(forward / thermal_voltage) + j02 * np.exp(forward / (2 * thermal_voltage))
        measured = ideal / (1.0 + resistance * ideal)
        axes[1, 0].semilogy(forward, measured, color=color, label=material)
    axes[1, 0].plot(forward, 1e-8 * np.exp(forward / thermal_voltage), "k:", label=r"$\eta=1$")
    axes[1, 0].plot(forward, 1e-5 * np.exp(forward / (2 * thermal_voltage)), "k--", label=r"$\eta=2$")
    axes[1, 0].set(xlabel="forward voltage (V)", ylabel="current (schematic)", title="Ideality and series-resistance regimes")
    axes[1, 0].legend(ncol=2)

    voltage_t = np.linspace(0.05, 0.82, 400)
    for temperature, color in ((250, BLUE), (300, GREEN), (400, ORANGE), (450, RED)):
        vt = K_B * temperature / Q
        gap = 1.17 - 4.73e-4 * temperature**2 / (temperature + 636)
        scale = temperature**3 * np.exp(-gap / vt)
        curve = scale * np.exp(voltage_t / vt)
        axes[1, 1].semilogy(voltage_t, curve / np.max(curve), color=color, label=f"{temperature} K")
    axes[1, 1].set(xlabel="forward voltage (V)", ylabel="normalized current", title="Temperature shifts forward characteristics")
    axes[1, 1].legend()
    for axis in axes.flat:
        _clean(axis)
    _save(figure, "junction-atlas-diode-characteristics.pdf")


def temperature_characteristics_atlas() -> None:
    """Generate forward and reverse silicon-diode temperature characteristics."""
    figure, axes = plt.subplots(1, 2, figsize=(9.2, 3.6), layout="constrained")
    forward_voltage = np.linspace(0.05, 0.82, 400)
    reverse_voltage = np.linspace(0.0, 12.0, 400)
    temperatures = ((250, BLUE), (300, GREEN), (400, ORANGE), (450, RED))
    intrinsic_scales = []
    for temperature, _ in temperatures:
        thermal_voltage = K_B * temperature / Q
        gap = 1.17 - 4.73e-4 * temperature**2 / (temperature + 636)
        intrinsic_scales.append(temperature**1.5 * np.exp(-gap / (2.0 * thermal_voltage)))
    reference = intrinsic_scales[1]
    for (temperature, color), intrinsic in zip(temperatures, intrinsic_scales, strict=True):
        thermal_voltage = K_B * temperature / Q
        gap = 1.17 - 4.73e-4 * temperature**2 / (temperature + 636)
        diffusion_scale = (intrinsic / reference) ** 2
        forward = diffusion_scale * np.exp((forward_voltage - 0.62) / thermal_voltage)
        generation = 2.0e-5 * (intrinsic / reference) * np.sqrt(0.8 + reverse_voltage)
        diffusion_reverse = 2.0e-7 * diffusion_scale * np.ones_like(reverse_voltage)
        reverse = generation + diffusion_reverse
        axes[0].semilogy(forward_voltage, forward, color=color, label=f"{temperature} K")
        axes[1].semilogy(reverse_voltage, reverse, color=color, label=f"{temperature} K")
    axes[0].set(xlabel="forward voltage (V)", ylabel="current magnitude (normalized)", title="Forward-bias temperature dependence")
    axes[1].set(xlabel="reverse-bias magnitude (V)", ylabel="current magnitude (normalized)", title="Generation-to-diffusion crossover")
    for axis in axes:
        axis.legend()
        _clean(axis)
    _save(figure, "junction-atlas-temperature-characteristics.pdf")


def small_signal_atlas() -> None:
    """Generate the diode small-signal conductance and capacitance circuit."""
    figure, axis = plt.subplots(figsize=(8.4, 3.0), layout="constrained")
    axis.axis("off")
    axis.plot([0.5, 7.8], [2.4, 2.4], color="black")
    axis.plot([0.5, 7.8], [0.4, 0.4], color="black")
    branches = ((2.0, r"$g_d$", "resistor"), (4.15, r"$C_j$", "capacitor"), (6.3, r"$C_d$", "capacitor"))
    for x, label, kind in branches:
        axis.plot([x, x], [2.4, 1.75], color="black")
        axis.plot([x, x], [1.05, 0.4], color="black")
        if kind == "capacitor":
            axis.plot([x - 0.35, x + 0.35], [1.62, 1.62], color=BLUE, linewidth=2)
            axis.plot([x - 0.35, x + 0.35], [1.18, 1.18], color=BLUE, linewidth=2)
        else:
            zigx = x + np.array([0, -0.25, 0.25, -0.25, 0.25, 0])
            zigy = np.linspace(1.75, 1.05, 6)
            axis.plot(zigx, zigy, color=RED)
        axis.text(x, 0.0, label, ha="center", fontsize=12)
    axis.text(0.35, 2.4, "+", fontsize=15, va="center")
    axis.text(0.35, 0.4, "−", fontsize=15, va="center")
    axis.text(4.15, 2.8, "linearized p–n-junction small-signal model", ha="center", fontweight="bold")
    axis.text(2.0, -0.45, "differential\nconductance", ha="center")
    axis.text(4.15, -0.45, "depletion\ncharge", ha="center")
    axis.text(6.3, -0.45, "stored minority\ncharge", ha="center")
    axis.set(xlim=(0, 8.2), ylim=(-0.7, 3.1))
    _save(figure, "junction-atlas-small-signal.pdf")


def switching_atlas() -> None:
    """Generate diode switching-circuit, recovery-waveform, and time panels."""
    figure, axes = plt.subplots(1, 3, figsize=(11.0, 3.5), layout="constrained")
    axes[0].axis("off")
    axes[0].plot([0.2, 1.1], [1.5, 1.5], color="black")
    axes[0].plot([1.1, 1.55], [1.5, 1.8], color="black", linewidth=2)
    axes[0].add_patch(Circle((1.7, 1.9), 0.06, color="black"))
    axes[0].add_patch(Circle((1.7, 1.1), 0.06, color="black"))
    axes[0].plot([1.75, 2.7], [1.9, 1.9], color=RED)
    axes[0].plot([1.75, 2.7], [1.1, 1.1], color=BLUE)
    axes[0].text(2.75, 1.9, r"$+V_F$")
    axes[0].text(2.75, 1.1, r"$-V_R$")
    axes[0].plot([0.2, 0.2], [1.15, 1.85], color="black")
    axes[0].plot([0.05, 0.35], [1.15, 1.15], color="black")
    axes[0].plot([0.05, 0.35], [1.85, 1.85], color="black")
    axes[0].text(0.2, 0.85, "diode", ha="center")
    axes[0].text(1.25, 2.25, "switch", ha="center")
    axes[0].set(xlim=(-0.1, 3.25), ylim=(0.55, 2.5), title="Forward-to-reverse switching circuit")

    time = np.linspace(-1.0, 5.0, 700)
    current = np.where(time < 0.0, 1.0, np.where(time < 2.0, -1.5, -1.5 * np.exp(-(time - 2.0) / 0.65)))
    axes[1].plot(time, current, color=BLUE)
    axes[1].axhline(0, color="black", linewidth=0.7)
    axes[1].axvspan(0, 2.0, color=ORANGE, alpha=0.18, label="storage")
    axes[1].axvspan(2.0, 4.2, color=GREEN, alpha=0.12, label="transition")
    axes[1].set(xlabel="time", ylabel="diode current", title="Reverse-recovery waveform")
    axes[1].legend()

    ratio = np.logspace(-2, 2, 300)
    reduced_time = np.log1p(ratio)
    axes[2].semilogx(ratio, reduced_time, color=GREEN)
    axes[2].set(xlabel=r"forward/reverse-current ratio $I_F/I_R$", ylabel=r"reduced removal time $t/\tau$", title="Charge-control switching trend")
    for axis in axes[1:]:
        _clean(axis)
    _save(figure, "junction-atlas-switching.pdf")


def breakdown_mechanisms_atlas() -> None:
    """Generate tunneling-band and avalanche-multiplication schematics."""
    figure, axes = plt.subplots(1, 3, figsize=(11.0, 3.5), layout="constrained")
    x = np.linspace(-1, 1, 300)
    conduction = 0.9 - 0.85 * (x + 1)
    valence = conduction - 0.9
    axes[0].plot(x, conduction, color=BLUE, label=r"$E_C$")
    axes[0].plot(x, valence, color=RED, label=r"$E_V$")
    axes[0].annotate("band-to-band\ntunneling", (0.18, -0.12), xytext=(-0.45, 0.35), arrowprops={"arrowstyle": "->", "color": GREEN}, color=GREEN)
    axes[0].set(xlabel="position", ylabel="energy", title="High-field tunneling")
    axes[0].set_yticks([])
    axes[0].legend()

    axes[1].plot(x, conduction, color=BLUE)
    axes[1].plot(x, valence, color=RED)
    starts = (-0.62, -0.15, 0.30)
    for index, start in enumerate(starts, 1):
        axes[1].annotate("", (start + 0.42, np.interp(start + 0.42, x, conduction) + 0.12), (start, np.interp(start, x, conduction) + 0.12), arrowprops={"arrowstyle": "->", "color": ORANGE, "lw": 1.6})
        axes[1].text(start + 0.16, np.interp(start, x, conduction) + 0.22, str(index), color=ORANGE)
    axes[1].annotate("new e–h pair", (0.23, -0.35), xytext=(-0.25, -0.75), arrowprops={"arrowstyle": "->", "color": GREEN}, color=GREEN)
    axes[1].set(xlabel="position", title="Impact-ionization cascade")
    axes[1].set_yticks([])

    distance = np.linspace(0, 1, 300)
    for integral, color in ((0.35, GREEN), (0.75, ORANGE), (0.98, RED)):
        multiplication = 1.0 / (1.0 - integral * distance)
        axes[2].plot(distance, multiplication, color=color, label=fr"$\int\alpha dx={integral}$")
    axes[2].set(xlabel=r"normalized depletion distance $x/W$", ylabel=r"multiplication $M(x)$", title="Current multiplication")
    axes[2].legend()
    for axis in axes:
        _clean(axis)
    _save(figure, "junction-atlas-breakdown-mechanisms.pdf")


def breakdown_design_atlas() -> None:
    """Generate critical-field, breakdown, diffused, and punch-through trends."""
    figure, axes = plt.subplots(2, 2, figsize=(9.2, 6.4), layout="constrained")
    doping_cm3 = np.logspace(14, 18.2, 50)
    doping_m3 = doping_cm3 * 1e6
    permittivity = 11.7 * EPSILON_0
    breakdown = np.array([avalanche_breakdown_voltage(n, 0.8, permittivity, 7e7, 1.23e8, maximum_voltage_v=2e5) for n in doping_m3])
    field = np.sqrt(2 * Q * doping_m3 * (breakdown + 0.8) / permittivity)
    axes[0, 0].loglog(doping_cm3, field / 1e2, color=BLUE, label="avalanche field")
    axes[0, 0].axhline(1e6, color=RED, linestyle="--", label="tunneling-scale field")
    axes[0, 0].set(xlabel=r"background doping (cm$^{-3}$)", ylabel="critical field (V cm$^{-1}$)", title="Critical-field competition")
    axes[0, 0].legend()

    axes[0, 1].loglog(doping_cm3, breakdown, color=BLUE, label="one-sided abrupt")
    axes[0, 1].loglog(doping_cm3, 1.7 * breakdown, color=GREEN, linestyle="--", label="graded trend")
    axes[0, 1].set(xlabel=r"doping or equivalent gradient scale", ylabel="breakdown voltage (V)", title="Avalanche design curves")
    axes[0, 1].legend()

    gradient_scale = np.logspace(-2, 2, 300)
    abrupt_limit = 25 * np.ones_like(gradient_scale)
    graded_limit = 70 / np.sqrt(gradient_scale)
    diffused = 1.0 / (1.0 / abrupt_limit + 1.0 / graded_limit)
    axes[1, 0].loglog(gradient_scale, abrupt_limit, label="abrupt limit")
    axes[1, 0].loglog(gradient_scale, graded_limit, label="graded limit")
    axes[1, 0].loglog(gradient_scale, diffused, linewidth=2, label="diffused junction")
    axes[1, 0].set(xlabel="profile parameter", ylabel="breakdown voltage", title="Diffused profile between limits")
    axes[1, 0].legend()

    doping_low = np.logspace(13, 17, 250)
    plane = 1.1e17 / doping_low
    for thickness, color in ((10, RED), (20, GREEN), (50, BLUE)):
        punch_limit = 18.0 * thickness
        axes[1, 1].loglog(doping_low, np.minimum(plane, punch_limit), color=color, label=fr"$W={thickness}\,\mu$m")
    axes[1, 1].set(xlabel=r"light-layer doping (cm$^{-3}$)", ylabel="supported voltage (V)", title="Punch-through voltage saturation")
    axes[1, 1].legend()
    for axis in axes.flat:
        _clean(axis)
    _save(figure, "junction-atlas-breakdown-design.pdf")


def curvature_atlas() -> None:
    """Generate planar-edge geometry and curvature-dependent breakdown panels."""
    figure, axes = plt.subplots(1, 2, figsize=(9.0, 3.6), layout="constrained")
    axes[0].add_patch(Rectangle((-1.6, 0.25), 3.2, 0.35, color="0.75", label="mask/oxide"))
    axes[0].add_patch(Rectangle((-0.75, 0.24), 1.5, 0.38, color="white"))
    theta = np.linspace(0, np.pi, 200)
    axes[0].plot(0.92 * np.cos(theta), -0.10 - 0.75 * np.sin(theta), color=RED, linewidth=2, label="diffused junction")
    axes[0].plot([-0.92, -0.92], [-0.1, 0.25], color=RED, linewidth=2)
    axes[0].plot([0.92, 0.92], [-0.1, 0.25], color=RED, linewidth=2)
    axes[0].annotate("field crowding", (0.87, -0.14), xytext=(1.15, -0.65), arrowprops={"arrowstyle": "->"})
    axes[0].set(xlim=(-1.8, 1.8), ylim=(-1.05, 0.8), xlabel="lateral coordinate", ylabel="depth", title="Planar diffusion creates curved edges")
    axes[0].legend(loc="lower left")
    axes[0].set_aspect("equal")

    normalized_radius = np.logspace(-1, 2, 300)
    axes[1].semilogx(normalized_radius, np.ones_like(normalized_radius), color="black", label="plane")
    axes[1].semilogx(normalized_radius, 1 - 0.38 / (1 + normalized_radius), color=BLUE, label="cylindrical")
    axes[1].semilogx(normalized_radius, 1 - 0.62 / (1 + normalized_radius), color=RED, label="spherical")
    axes[1].set(xlabel=r"junction radius / plane depletion width", ylabel=r"$V_B/V_{B,\mathrm{plane}}$", title="Curvature lowers breakdown voltage")
    axes[1].legend()
    _clean(axes[1])
    _save(figure, "junction-atlas-curvature.pdf")


def heterojunction_atlas() -> None:
    """Generate isolated-material and equilibrium heterojunction band diagrams."""
    figure, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), layout="constrained")
    for axis in axes:
        axis.set_yticks([])
        _clean(axis)
    axes[0].axhline(1.45, color="black", label=r"$E_{vac}$")
    axes[0].plot([-1.0, -0.1], [0.0, 0.0], color=BLUE)
    axes[0].plot([-1.0, -0.1], [-1.42, -1.42], color=RED)
    axes[0].plot([0.1, 1.0], [0.57, 0.57], color=BLUE)
    axes[0].plot([0.1, 1.0], [-1.59, -1.59], color=RED)
    axes[0].plot([-1.0, -0.1], [-0.20, -0.20], color=GRAY, linestyle="--")
    axes[0].plot([0.1, 1.0], [-1.39, -1.39], color=GRAY, linestyle="--")
    axes[0].annotate(r"$\chi_1$", (-0.7, 0.0), (-0.7, 1.25), arrowprops={"arrowstyle": "<->"}, ha="center")
    axes[0].annotate(r"$\chi_2$", (0.7, 0.57), (0.7, 1.25), arrowprops={"arrowstyle": "<->"}, ha="center")
    axes[0].text(-0.55, -0.25, "material 1", ha="center")
    axes[0].text(0.55, -0.95, "material 2", ha="center")
    axes[0].set(xlim=(-1.1, 1.1), ylim=(-2.45, 1.65), xlabel="isolated semiconductors", title="Before contact: affinity, gap, and work function")

    x = np.linspace(-1, 1, 400)
    bend = 0.85 * (1 + np.tanh(3 * x)) / 2
    ec = bend + np.where(x >= 0, 0.57, 0.0)
    ev = np.where(x < 0, ec - 1.42, ec - 2.16)
    axes[1].plot(x, ec, color=BLUE, label=r"$E_C$")
    axes[1].plot(x, ev, color=RED, label=r"$E_V$")
    axes[1].axhline(0, color=GRAY, linestyle="--", label=r"$E_F$")
    axes[1].axvline(0, color="0.55", linestyle=":")
    axes[1].annotate(r"$\Delta E_C$", (0.03, ec[len(x)//2]), (0.35, 1.45), arrowprops={"arrowstyle": "->"})
    axes[1].annotate(r"$\Delta E_V$", (0.03, ev[len(x)//2]), (-0.45, -1.25), arrowprops={"arrowstyle": "->"})
    axes[1].set(xlim=(-1, 1), ylim=(-1.7, 1.75), xlabel="position", title="After contact: flat Fermi level and band offsets")
    axes[1].legend()
    _save(figure, "junction-atlas-heterojunction.pdf")


def generate_all_junction_atlas_figures() -> None:
    """Generate every original schematic plate for Sze--Lee Chapter 3 topics."""
    junction_formation_atlas()
    doping_electrostatics_atlas()
    graded_junction_atlas()
    capacitance_construction_atlas()
    injection_transport_atlas()
    diode_characteristics_atlas()
    temperature_characteristics_atlas()
    small_signal_atlas()
    switching_atlas()
    breakdown_mechanisms_atlas()
    breakdown_design_atlas()
    curvature_atlas()
    heterojunction_atlas()
