"""Reproducible equilibrium-semiconductor publication figures.

This module contains the fixed scientific recipes used by the intrinsic and
extrinsic semiconductor chapters.  The figures distinguish available states,
occupation probabilities, carrier populations, dopant charge states, and
charge-neutrality closure.  Their models are intentionally compact enough to
inspect: parabolic bands, explicit Fermi--Dirac or Boltzmann occupations,
standard donor and acceptor degeneracy factors, and scalar neutrality solves.

Each public recipe has no arguments because its numerical values are part of
the manuscript's declared example.  It writes one timestamp-free PDF beneath
``figures/`` and closes the Matplotlib object.  The corresponding student
notebook first evaluates the governing relation directly and imports the
recipe only at the end.  These routines support teaching and reproducibility;
they are not calibrated material simulators or publication-grade device
software.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

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


def save(fig: plt.Figure, name: str) -> None:
    """Write one deterministic PDF and close its Matplotlib figure.

    Parameters
    ----------
    fig:
        Completed Matplotlib figure.  Ownership transfers to this function,
        which closes the object after writing it.
    name:
        PDF basename relative to :data:`FIGURE_DIR`.

    Returns
    -------
    None
        The generated PDF is the function's external result.

    Notes
    -----
    Creation and modification timestamps are removed from the PDF metadata so
    repeated executions can be compared byte for byte.
    """
    metadata = {
        "Creator": "SEMPHY01 deterministic figure generator",
        "Producer": "Matplotlib",
        "CreationDate": None,
        "ModDate": None,
    }
    fig.savefig(FIGURE_DIR / name, metadata=metadata)
    plt.close(fig)


def band_formation_schematic() -> None:
    """Show how repeated atomic levels become families of crystal states.

    Returns
    -------
    None
        Writes ``intrinsic-band-formation.pdf``.

    Notes
    -----
    The vertical axis is schematic.  Two isolated one-electron levels split as
    equivalent atoms interact and are finally represented as valence and
    conduction bands separated by a forbidden gap.  The drawing communicates
    spectral multiplicity rather than a material-specific dispersion.
    """
    fig, axis = plt.subplots(figsize=(7.8, 4.2), layout="constrained")
    axis.set_xlim(-0.5, 2.8)
    axis.set_ylim(-0.4, 3.8)
    axis.set_xticks([0.0, 1.1, 2.25], ["isolated atom", "interacting atoms", "crystal"])
    axis.set_ylabel("one-electron energy (schematic)")
    axis.set_yticks([])

    atomic_levels = [0.7, 2.6]
    for energy_level in atomic_levels:
        axis.hlines(energy_level, -0.18, 0.18, color="#123B66", lw=2.2)

    lower_split = np.linspace(0.45, 0.95, 8)
    upper_split = np.linspace(2.35, 2.85, 8)
    for energy_level in np.concatenate((lower_split, upper_split)):
        axis.hlines(energy_level, 0.80, 1.40, color="#4C78A8", lw=0.8)

    axis.fill_between(
        [1.85, 2.65], [0.35, 0.35], [1.10, 1.10], color="#4C78A8", alpha=0.35
    )
    axis.fill_between(
        [1.85, 2.65], [2.20, 2.20], [3.00, 3.00], color="#F28E2B", alpha=0.35
    )
    axis.hlines([0.35, 1.10], 1.85, 2.65, color="#123B66", lw=1.2)
    axis.hlines([2.20, 3.00], 1.85, 2.65, color="#A44A17", lw=1.2)
    axis.annotate(
        r"forbidden gap $E_g$",
        xy=(2.25, 1.65),
        ha="center",
        va="center",
        color="#A44A17",
    )
    axis.text(2.25, 0.72, "valence band", ha="center", va="center")
    axis.text(2.25, 2.60, "conduction band", ha="center", va="center")
    axis.annotate(
        "interatomic coupling\nand Pauli exclusion",
        xy=(1.70, 3.42),
        xytext=(0.45, 3.42),
        arrowprops={"arrowstyle": "->", "color": "0.3"},
        ha="center",
        va="center",
    )
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.grid(False)
    axis.set_title("From atomic levels to allowed energy bands")
    save(fig, "intrinsic-band-formation.pdf")


def kspace_state_counting() -> None:
    """Connect a finite reciprocal-space grid to a continuum band-edge DOS.

    Returns
    -------
    None
        Writes ``intrinsic-kspace-state-counting.pdf``.

    Notes
    -----
    The first panel counts grid points inside an elliptical constant-energy
    contour, representing a finite periodic cell.  The second panel shows the
    three-dimensional parabolic laws ``N(E) proportional to E^(3/2)`` and
    ``g(E) proportional to sqrt(E)``.  Axes are normalized because the purpose
    is to expose counting and differentiation rather than assign a material
    prefactor.
    """
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.6), layout="constrained")

    q1 = np.linspace(-2.5, 2.5, 19)
    q2 = np.linspace(-2.0, 2.0, 17)
    mesh_q1, mesh_q2 = np.meshgrid(q1, q2)
    scaled_energy = (mesh_q1 / 2.1) ** 2 + (mesh_q2 / 1.15) ** 2
    inside = scaled_energy <= 1.0
    axes[0].scatter(
        mesh_q1[~inside], mesh_q2[~inside], s=7, color="0.82", label="higher energy"
    )
    axes[0].scatter(
        mesh_q1[inside], mesh_q2[inside], s=11, color="#123B66", label=r"$E\leq E_C+\varepsilon$"
    )
    angle = np.linspace(0.0, 2.0 * np.pi, 400)
    axes[0].plot(2.1 * np.cos(angle), 1.15 * np.sin(angle), color="#A44A17", lw=1.6)
    axes[0].axhline(0.0, color="0.5", lw=0.6)
    axes[0].axvline(0.0, color="0.5", lw=0.6)
    axes[0].set_aspect("equal")
    axes[0].set_xlabel(r"$q_1$")
    axes[0].set_ylabel(r"$q_2$")
    axes[0].set_title("finite periodic cell")
    axes[0].legend(frameon=False, loc="upper right")

    epsilon = np.linspace(0.0, 1.0, 400)
    integrated = epsilon**1.5
    density = np.sqrt(epsilon)
    axes[1].plot(epsilon, integrated, color="#123B66", lw=2, label=r"$\mathcal{N}(\varepsilon)\propto\varepsilon^{3/2}$")
    axes[1].plot(epsilon, density, color="#A44A17", lw=2, label=r"$g(\varepsilon)\propto\sqrt{\varepsilon}$")
    axes[1].set_xlabel(r"edge energy $\varepsilon$")
    axes[1].set_ylabel("normalized count")
    axes[1].set_title("continuum state count")
    axes[1].legend(frameon=False)
    axes[1].grid(alpha=0.2)

    fig.suptitle("From discrete wave vectors to a band-edge density of states")
    save(fig, "intrinsic-kspace-state-counting.pdf")


def fermi_boltzmann_comparison() -> None:
    """Compare exact Fermi--Dirac occupation with its Boltzmann limit.

    Returns
    -------
    None
        Writes ``intrinsic-fermi-boltzmann.pdf``.

    Notes
    -----
    Energy is measured by the dimensionless variable ``(E-mu)/(k_B T)``.  The
    left panel overlays the exact occupation and its exponential tail.  The
    right panel plots the relative error, which equals ``exp[-(E-mu)/(k_B T)]``
    for this comparison and reaches approximately five percent at three
    thermal energies above the chemical potential.
    """
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.5), layout="constrained")
    reduced_energy = np.linspace(-6.0, 8.0, 800)
    exact = 1.0 / (1.0 + np.exp(reduced_energy))
    boltzmann = np.exp(-reduced_energy)

    axes[0].plot(reduced_energy, exact, color="#123B66", lw=2, label="Fermi--Dirac")
    positive = reduced_energy >= 0.0
    axes[0].plot(
        reduced_energy[positive],
        boltzmann[positive],
        color="#A44A17",
        ls="--",
        lw=1.8,
        label="Boltzmann tail",
    )
    axes[0].axvline(3.0, color="0.45", ls=":", lw=1)
    axes[0].set_xlabel(r"$(E-\mu)/(k_{\mathrm{B}}T)$")
    axes[0].set_ylabel("occupation probability")
    axes[0].set_ylim(-0.03, 1.05)
    axes[0].set_title("occupation")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.2)

    tail_energy = np.linspace(0.0, 8.0, 500)
    relative_error = np.exp(-tail_energy)
    axes[1].semilogy(tail_energy, relative_error, color="#2E7D32", lw=2)
    axes[1].axvline(3.0, color="0.45", ls=":", lw=1)
    axes[1].axhline(np.exp(-3.0), color="0.45", ls=":", lw=1)
    axes[1].annotate(
        r"at $3k_{\mathrm{B}}T$: $e^{-3}\approx5\%$",

        xy=(3.0, np.exp(-3.0)),
        xytext=(4.2, 0.16),
        arrowprops={"arrowstyle": "->", "color": "0.35"},
    )
    axes[1].set_xlabel(r"$(E-\mu)/(k_{\mathrm{B}}T)$")
    axes[1].set_ylabel("relative error of Boltzmann approximation")
    axes[1].set_title("approximation error")
    axes[1].grid(alpha=0.2, which="both")
    save(fig, "intrinsic-fermi-boltzmann.pdf")


def states_occupation_population() -> None:
    """Separate available states, occupation, and carrier spectral density.

    Returns
    -------
    None
        Writes ``intrinsic-states-occupation.pdf``.

    Notes
    -----
    Square-root valence and conduction densities of states form an asymmetric
    illustrative semiconductor.  A bisection solve selects the chemical
    potential for intrinsic neutrality.  The panels then display ``g(E)``,
    ``f(E)`` and ``1-f(E)``, and the products ``g_c f`` and ``g_v(1-f)`` whose
    integrals are the electron and hole populations.  An assertion verifies
    equality of those integrals to relative tolerance ``1e-12``.
    """
    ev = -0.55
    ec = 0.55
    energy = np.linspace(-1.05, 1.05, 1600)
    thermal_energy = 0.10

    conduction_dos = np.sqrt(np.clip(energy - ec, 0.0, None))
    valence_dos = 1.15 * np.sqrt(np.clip(ev - energy, 0.0, None))

    # Solve intrinsic neutrality for the deliberately asymmetric DOS model.
    lower_mu = ev
    upper_mu = ec
    for _ in range(100):
        trial_mu = 0.5 * (lower_mu + upper_mu)
        trial_occupation = 1.0 / (
            1.0 + np.exp((energy - trial_mu) / thermal_energy)
        )
        electron_count = np.trapezoid(
            conduction_dos * trial_occupation, energy
        )
        hole_count = np.trapezoid(
            valence_dos * (1.0 - trial_occupation), energy
        )
        if electron_count > hole_count:
            upper_mu = trial_mu
        else:
            lower_mu = trial_mu
    chemical_potential = 0.5 * (lower_mu + upper_mu)
    occupation = 1.0 / (
        1.0 + np.exp((energy - chemical_potential) / thermal_energy)
    )
    electrons = conduction_dos * occupation
    holes = valence_dos * (1.0 - occupation)
    assert np.isclose(
        np.trapezoid(electrons, energy),
        np.trapezoid(holes, energy),
        rtol=1.0e-12,
    )

    fig, axes = plt.subplots(
        1, 3, figsize=(9.2, 3.4), sharey=True, layout="constrained"
    )

    axes[0].plot(conduction_dos, energy, color="#1f77b4", label=r"$g_c(E)$")
    axes[0].plot(valence_dos, energy, color="#d62728", label=r"$g_v(E)$")
    axes[0].set_xlabel("density of states (arb. units)")
    axes[0].set_ylabel("energy")
    axes[0].set_title("available states")
    axes[0].legend(loc="center right", frameon=False)

    axes[1].plot(occupation, energy, color="#2ca02c", label=r"$f(E)$")
    axes[1].plot(1.0 - occupation, energy, color="#9467bd", label=r"$1-f(E)$")
    axes[1].axvline(0.5, color="0.75", lw=0.8)
    axes[1].set_xlim(-0.02, 1.02)
    axes[1].set_xlabel("probability")
    axes[1].set_title("occupation")
    axes[1].legend(loc="center right", frameon=False)

    axes[2].plot(electrons, energy, color="#1f77b4", label=r"$g_cf$")
    axes[2].fill_betweenx(energy, 0.0, electrons, color="#1f77b4", alpha=0.22)
    axes[2].plot(holes, energy, color="#d62728", label=r"$g_v(1-f)$")
    axes[2].fill_betweenx(energy, 0.0, holes, color="#d62728", alpha=0.22)
    axes[2].set_xlabel("carrier spectral density (arb. units)")
    axes[2].set_title("populations")
    axes[2].legend(loc="center right", frameon=False)

    for axis in axes:
        axis.axhspan(ev, ec, color="0.93", zorder=-10)
        axis.axhline(ec, color="0.35", ls="--", lw=0.8)
        axis.axhline(ev, color="0.35", ls="--", lw=0.8)
        axis.axhline(chemical_potential, color="black", ls=":", lw=0.9)
        axis.grid(alpha=0.15)

    y_axis_transform = axes[2].get_yaxis_transform()
    axes[2].text(
        0.98,
        ec + 0.03,
        r"$E_C$",
        ha="right",
        va="bottom",
        transform=y_axis_transform,
    )
    axes[2].text(
        0.98,
        ev - 0.03,
        r"$E_V$",
        ha="right",
        va="top",
        transform=y_axis_transform,
    )
    axes[2].text(
        0.98,
        chemical_potential + 0.03,
        r"$\mu_i$",
        ha="right",
        transform=y_axis_transform,
    )
    fig.suptitle("Intrinsic semiconductor: states, occupation, and carriers")
    save(fig, "intrinsic-states-occupation.pdf")


def intrinsic_temperature_dependence() -> None:
    """Evaluate an illustrative intrinsic-silicon temperature model.

    Returns
    -------
    None
        Writes ``intrinsic-temperature-dependence.pdf``.

    Notes
    -----
    The recipe combines a Varshni-type gap, ``T^(3/2)`` effective densities of
    states, and the nondegenerate intrinsic relation
    ``n_i = sqrt(N_c N_v) exp[-E_g/(2 k_B T)]``.  The second panel displays the
    mass-asymmetry shift of the intrinsic chemical potential away from midgap.
    Numerical parameters are the fixed illustrative values stated in the
    manuscript caption, not a fitted data product.
    """
    kb_ev = 8.617333262145e-5
    temperature = np.linspace(200.0, 600.0, 401)

    # Illustrative silicon parameters following the conventions discussed in
    # the chapter. They are displayed in the figure caption, not hidden.
    eg0_ev = 1.169
    alpha_ev_per_k = 4.9e-4
    beta_k = 655.0
    eg_ev = eg0_ev - alpha_ev_per_k * temperature**2 / (temperature + beta_k)

    nc_300 = 2.86e19
    nv_300 = 2.66e19
    nc = nc_300 * (temperature / 300.0) ** 1.5
    nv = nv_300 * (temperature / 300.0) ** 1.5
    intrinsic = np.sqrt(nc * nv) * np.exp(
        -eg_ev / (2.0 * kb_ev * temperature)
    )
    mu_shift = 0.5 * kb_ev * temperature * np.log(nv / nc)

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.5), layout="constrained")
    axes[0].semilogy(temperature, intrinsic, color="#123B66", lw=2)
    axes[0].axvline(300.0, color="0.4", ls="--", lw=0.8)
    axes[0].set_xlabel("temperature (K)")
    axes[0].set_ylabel(r"$n_i$ (cm$^{-3}$)")
    axes[0].set_title("intrinsic concentration")
    axes[0].grid(alpha=0.2, which="both")

    axes[1].plot(temperature, 1.0e3 * mu_shift, color="#A44A17", lw=2)
    axes[1].axhline(0.0, color="0.4", ls="--", lw=0.8)
    axes[1].set_xlabel("temperature (K)")
    axes[1].set_ylabel(r"$\mu_i-E_{\mathrm{mid}}$ (meV)")
    axes[1].set_title("intrinsic chemical-potential shift")
    axes[1].grid(alpha=0.2)

    fig.suptitle("Illustrative parabolic-band silicon model")
    save(fig, "intrinsic-temperature-dependence.pdf")


def extrinsic_levels_and_charge_states() -> None:
    """Draw localized donor and acceptor levels with their charge bookkeeping.

    Returns
    -------
    None
        Writes ``extrinsic-dopant-levels.pdf``.

    Notes
    -----
    A shallow donor level lies below the conduction-band edge and leaves a
    positive ion when its electron enters the conduction band.  A shallow
    acceptor level lies above the valence-band edge and becomes negative when
    it receives a valence electron, leaving a mobile hole.  Energies are
    schematic fractions of the displayed gap.
    """
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6), sharey=True, layout="constrained")
    ev = 0.0
    ec = 1.0
    ed = 0.84
    ea = 0.16

    for axis, title in zip(axes, ("donor", "acceptor"), strict=True):
        axis.set_xlim(0.0, 1.0)
        axis.set_ylim(-0.12, 1.12)
        axis.axhspan(ev, ec, color="0.95")
        axis.hlines(ec, 0.08, 0.92, color="#123B66", lw=2)
        axis.hlines(ev, 0.08, 0.92, color="#123B66", lw=2)
        axis.text(0.94, ec, r"$E_C$", va="center", ha="right")
        axis.text(0.94, ev, r"$E_V$", va="center", ha="right")
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_title(title)
        axis.spines[:].set_visible(False)

    axes[0].hlines(ed, 0.18, 0.72, color="#A44A17", lw=2)
    axes[0].text(0.75, ed, r"$E_D$", va="center")
    axes[0].annotate(
        r"electron to $E_C$",
        xy=(0.44, ec - 0.01),
        xytext=(0.44, ed + 0.02),
        arrowprops={"arrowstyle": "->", "color": "#A44A17"},
        ha="center",
        va="bottom",
    )
    axes[0].text(0.44, 0.47, r"$D^0\rightarrow D^+ + e^-$", ha="center")

    axes[1].hlines(ea, 0.18, 0.72, color="#2E7D32", lw=2)
    axes[1].text(0.75, ea, r"$E_A$", va="center")
    axes[1].annotate(
        r"valence electron to $E_A$",
        xy=(0.44, ea + 0.01),
        xytext=(0.44, ev + 0.01),
        arrowprops={"arrowstyle": "->", "color": "#2E7D32"},
        ha="center",
        va="bottom",
    )
    axes[1].text(0.44, 0.58, r"$A^0\rightarrow A^- + h^+$", ha="center")
    fig.suptitle("Localized dopant levels and mobile carriers")
    save(fig, "extrinsic-dopant-levels.pdf")


def extrinsic_ionization_fractions() -> None:
    """Plot donor and acceptor ionization with explicit degeneracy factors.

    Returns
    -------
    None
        Writes ``extrinsic-ionization-fractions.pdf``.

    Notes
    -----
    The donor curve uses ``N_D^+/N_D = 1/[1+2 exp((mu-E_D)/k_B T)]`` and the
    acceptor curve uses ``N_A^-/N_A = 1/[1+4 exp((E_A-mu)/k_B T)]`` in the
    Sze--Ng convention adopted by the manuscript.  Dashed guides identify the
    degeneracy-shifted half-ionization points.
    """
    reduced_mu = np.linspace(-8.0, 8.0, 800)
    donor_fraction = 1.0 / (1.0 + 2.0 * np.exp(reduced_mu))
    acceptor_fraction = 1.0 / (1.0 + 4.0 * np.exp(-reduced_mu))

    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.3), layout="constrained")
    axes[0].plot(reduced_mu, donor_fraction, color="#A44A17", lw=2)
    axes[0].set_xlabel(r"$(\mu-E_D)/(k_{\mathrm{B}}T)$")
    axes[0].set_ylabel(r"ionized fraction $N_D^+/N_D$")
    axes[0].set_title("donor ionization")

    axes[1].plot(reduced_mu, acceptor_fraction, color="#2E7D32", lw=2)
    axes[1].set_xlabel(r"$(\mu-E_A)/(k_{\mathrm{B}}T)$")
    axes[1].set_ylabel(r"ionized fraction $N_A^-/N_A$")
    axes[1].set_title("acceptor ionization")

    for axis in axes:
        axis.set_ylim(-0.03, 1.03)
        axis.axhline(0.5, color="0.55", ls=":", lw=0.8)
        axis.axvline(0.0, color="0.55", ls=":", lw=0.8)
        axis.grid(alpha=0.2)
    axes[0].axvline(-np.log(2.0), color="#A44A17", ls="--", lw=0.9)
    axes[1].axvline(np.log(4.0), color="#2E7D32", ls="--", lw=0.9)
    save(fig, "extrinsic-ionization-fractions.pdf")


def extrinsic_temperature_regimes() -> None:
    """Solve donor-doped charge neutrality across temperature regimes.

    Returns
    -------
    None
        Writes ``extrinsic-temperature-regimes.pdf``.

    Notes
    -----
    At each temperature, bisection determines the chemical potential satisfying
    ``n-p-N_D^+=0``.  The model includes incomplete donor ionization, a
    temperature-dependent gap, and parabolic effective densities of states.
    The plotted solution exposes freeze-out, the extrinsic plateau, and the
    intrinsic regime.  A normalized residual assertion verifies neutrality to
    better than ``1e-10`` before the figure is written.
    """
    kb_ev = 8.617333262145e-5
    temperature = np.linspace(30.0, 1000.0, 800)
    donor_concentration = 1.0e15
    donor_degeneracy = 2.0
    donor_binding_ev = 0.045

    eg0_ev = 1.169
    alpha_ev_per_k = 4.9e-4
    beta_k = 655.0
    gap = eg0_ev - alpha_ev_per_k * temperature**2 / (temperature + beta_k)
    ev = np.zeros_like(temperature)
    ec = gap
    ed = ec - donor_binding_ev
    nc = 2.86e19 * (temperature / 300.0) ** 1.5
    nv = 2.66e19 * (temperature / 300.0) ** 1.5

    electron = np.empty_like(temperature)
    hole = np.empty_like(temperature)
    ionized = np.empty_like(temperature)
    for index, temp in enumerate(temperature):
        thermal = kb_ev * temp
        lower_mu = ev[index] - 0.25
        upper_mu = ec[index] + 0.25
        for _ in range(160):
            mu = 0.5 * (lower_mu + upper_mu)
            n_value = nc[index] * np.exp((mu - ec[index]) / thermal)
            p_value = nv[index] * np.exp((ev[index] - mu) / thermal)
            donor_plus = donor_concentration / (
                1.0
                + donor_degeneracy * np.exp((mu - ed[index]) / thermal)
            )
            residual = n_value - p_value - donor_plus
            if residual > 0.0:
                upper_mu = mu
            else:
                lower_mu = mu
        mu = 0.5 * (lower_mu + upper_mu)
        electron[index] = nc[index] * np.exp((mu - ec[index]) / thermal)
        hole[index] = nv[index] * np.exp((ev[index] - mu) / thermal)
        ionized[index] = donor_concentration / (
            1.0 + donor_degeneracy * np.exp((mu - ed[index]) / thermal)
        )

    assert np.max(np.abs(electron - hole - ionized) / electron) < 1.0e-10
    intrinsic = np.sqrt(nc * nv) * np.exp(-gap / (2.0 * kb_ev * temperature))

    fig, axis = plt.subplots(figsize=(7.4, 3.8), layout="constrained")
    inverse_temperature = 1000.0 / temperature
    axis.semilogy(inverse_temperature, electron, color="#123B66", lw=2, label=r"$n$")
    axis.semilogy(
        inverse_temperature,
        np.where(hole >= 1.0e7, hole, np.nan),
        color="#A44A17",
        lw=1.8,
        label=r"$p$",
    )
    axis.semilogy(
        inverse_temperature,
        np.where(intrinsic >= 1.0e7, intrinsic, np.nan),
        color="#7B3F98",
        lw=1.4,
        ls="-.",
        label=r"$n_i$",
    )
    axis.semilogy(
        inverse_temperature,
        ionized,
        color="#2E7D32",
        lw=1.5,
        ls="--",
        label=r"$N_D^+$",
    )
    axis.axhline(donor_concentration, color="0.45", ls=":", lw=0.9, label=r"$N_D$")
    axis.text(27.0, 3.0e18, "freeze-out", ha="center", color="#A44A17")
    axis.text(8.0, 3.0e18, "extrinsic", ha="center", color="#2E7D32")
    axis.text(1.55, 3.0e18, "intrinsic", ha="center", color="#7B3F98")
    axis.invert_xaxis()
    axis.set_ylim(1.0e7, 1.0e19)
    axis.set_xlabel(r"$1000/T$ (K$^{-1}$)")
    axis.set_ylabel(r"concentration (cm$^{-3}$)")
    axis.set_title("Illustrative donor-doped temperature regimes")
    axis.legend(frameon=False, ncol=5, loc="lower left")
    axis.grid(alpha=0.2, which="both")
    save(fig, "extrinsic-temperature-regimes.pdf")


def extrinsic_compensation() -> None:
    """Solve the complete-ionization compensation model analytically.

    Returns
    -------
    None
        Writes ``extrinsic-compensation.pdf``.

    Notes
    -----
    Charge neutrality ``n-p=N_D-N_A`` and mass action ``np=n_i^2`` reduce to a
    quadratic equation.  Its positive roots give electron and hole populations
    continuously across acceptor-dominated, intrinsic, and donor-dominated
    conditions.  A symmetric logarithmic horizontal scale preserves the
    compensated neighborhood near zero net doping.
    """
    intrinsic = 1.0e10
    net_doping = np.linspace(-1.0e17, 1.0e17, 2001)
    root = np.sqrt(net_doping**2 + 4.0 * intrinsic**2)
    electron = 0.5 * (net_doping + root)
    hole = 0.5 * (-net_doping + root)

    fig, axis = plt.subplots(figsize=(7.4, 3.7), layout="constrained")
    axis.semilogy(net_doping, electron, color="#123B66", lw=2, label=r"$n$")
    axis.semilogy(net_doping, hole, color="#A44A17", lw=2, label=r"$p$")
    axis.axvline(0.0, color="0.45", ls=":", lw=0.9)
    axis.axhline(intrinsic, color="#2E7D32", ls="--", lw=1, label=r"$n_i$")
    axis.set_xscale("symlog", linthresh=1.0e10)
    axis.set_xticks([-1.0e17, -1.0e14, -1.0e11, 0.0, 1.0e11, 1.0e14, 1.0e17])
    axis.set_xlabel(r"net ionized doping $N_D-N_A$ (cm$^{-3}$)")
    axis.set_ylabel(r"carrier concentration (cm$^{-3}$)")
    axis.set_title("Complete-ionization, nondegenerate compensation model")
    axis.legend(frameon=False)
    axis.grid(alpha=0.2, which="both")
    save(fig, "extrinsic-compensation.pdf")


def main() -> None:
    """Regenerate every intrinsic and extrinsic figure in this module.

    Returns
    -------
    None
        Writes nine PDFs into :data:`FIGURE_DIR`.

    Notes
    -----
    The teaching workflow normally executes one recipe from each figure's
    notebook.  This aggregate entry point is retained for batch verification
    and deterministic hash checks.
    """
    band_formation_schematic()
    kspace_state_counting()
    fermi_boltzmann_comparison()
    states_occupation_population()
    intrinsic_temperature_dependence()
    extrinsic_levels_and_charge_states()
    extrinsic_ionization_fractions()
    extrinsic_temperature_regimes()
    extrinsic_compensation()


if __name__ == "__main__":
    main()
