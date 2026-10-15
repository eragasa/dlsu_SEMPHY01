"""Publication figure recipes for spectral density-of-states concepts.

The functions in this module translate the reusable numerical relations in
:mod:`semphy.density_of_states` into the three figures used by Chapter 3.  They
cover the inverse-slope origin of the one-dimensional density of states, the
dimensional power law of a parabolic band, and convergence with reciprocal-
space sampling and energy broadening.

Every recipe fixes its physical parameters, sampling grid, labels, and PDF
metadata.  A recipe therefore represents a published calculation rather than
a general plotting interface.  Students encounter the explicit NumPy
construction first in the corresponding one-figure notebook and call these
functions only in its final reusable-implementation section.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from .density_of_states import (
    ELECTRON_MASS_KG,
    parabolic_density_of_states,
    sampled_band_density_of_states,
    tight_binding_band_1d,
    tight_binding_dos_1d_analytic,
    uniform_bz_grid_1d,
)

FIGURE_DIR = Path(__file__).resolve().parents[2]

mpl.rcParams.update(
    {
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "legend.fontsize": 8,
        "figure.dpi": 120,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
        "pdf.fonttype": 42,
    }
)

METADATA = {
    "Creator": "SEMPHY01 deterministic DOS figure generator",
    "Producer": "Matplotlib",
    "CreationDate": None,
    "ModDate": None,
}


def save(figure: plt.Figure, name: str) -> None:
    """Write a timestamp-free PDF and release its Matplotlib resources.

    Parameters
    ----------
    figure:
        Completed Matplotlib figure to be written and then closed.
    name:
        Basename of the PDF under :data:`FIGURE_DIR`.

    Returns
    -------
    None
        The function communicates through the generated PDF file.

    Notes
    -----
    Timestamp fields are deliberately ``None``.  This makes repeated notebook
    executions byte reproducible when the numerical and plotting dependencies
    are unchanged.
    """
    figure.savefig(FIGURE_DIR / name, metadata=METADATA)
    plt.close(figure)


def band_to_dos() -> None:
    """Relate a cosine Bloch band to its analytic one-dimensional DOS.

    The left panel evaluates ``E(k) = -2t cos(ka)`` over the first Brillouin
    zone.  The right panel evaluates the corresponding spin-degenerate DOS per
    primitive cell.  Shared energy axes show that energies with a small band
    slope accumulate more states and develop van Hove singularities at the
    band edges.

    Returns
    -------
    None
        Writes ``dos-band-to-density.pdf`` into :data:`FIGURE_DIR`.
    """
    lattice_constant = 0.5e-9
    reduced_k = np.linspace(-1.0, 1.0, 2001)
    wave_vector = reduced_k * np.pi / lattice_constant
    energy = tight_binding_band_1d(
        wave_vector,
        lattice_constant_m=lattice_constant,
        hopping_ev=1.0,
    )
    energy_grid = np.linspace(-1.995, 1.995, 2000)
    dos = tight_binding_dos_1d_analytic(
        energy_grid,
        hopping_ev=1.0,
        degeneracy=2.0,
    )

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(8.4, 3.5),
        sharey=True,
        gridspec_kw={"width_ratios": [1.2, 0.8]},
        layout="constrained",
    )
    axes[0].plot(reduced_k, energy, color="#123B66", lw=2)
    axes[0].set_xlabel(r"reduced wave vector $ka/\pi$")
    axes[0].set_ylabel(r"energy $E$ (eV)")
    axes[0].set_title(r"band $E(k)=-2t\cos(ka)$")
    axes[0].axvline(0.0, color="0.6", ls=":", lw=0.8)

    axes[1].plot(dos, energy_grid, color="#A44A17", lw=2)
    axes[1].set_xlabel(r"DOS (states eV$^{-1}$ cell$^{-1}$)")
    axes[1].set_title("states accumulated at each energy")
    axes[1].set_xlim(0.0, 3.0)

    for selected_energy in (-1.5, 0.0, 1.5):
        for axis in axes:
            axis.axhline(selected_energy, color="0.7", ls="--", lw=0.7)
    axes[1].annotate(
        "small band slope\n$\\Rightarrow$ large DOS",
        xy=(2.25, 1.85),
        xytext=(1.0, 0.9),
        arrowprops={"arrowstyle": "->", "color": "0.35"},
        ha="center",
    )
    save(fig, "dos-band-to-density.pdf")


def dimensionality_comparison() -> None:
    """Compare parabolic-band DOS power laws in one, two, and three dimensions.

    Each curve is calculated from the continuum state-counting expression
    ``g_d(E) proportional to (E-E0)^(d/2-1)`` using the same effective mass and
    spin degeneracy.  Normalization at 0.25 eV removes the dimension-dependent
    units while retaining the distinct edge shapes: divergent in 1D, constant
    in 2D, and square-root in 3D.

    Returns
    -------
    None
        Writes ``dos-dimensionality.pdf`` into :data:`FIGURE_DIR`.
    """
    excess = np.linspace(0.002, 0.5, 1000)
    reference_index = np.argmin(np.abs(excess - 0.25))
    colors = {1: "#A44A17", 2: "#2E7D32", 3: "#123B66"}

    fig, axis = plt.subplots(figsize=(6.8, 3.7), layout="constrained")
    for dimension in (1, 2, 3):
        dos = parabolic_density_of_states(
            excess,
            edge_ev=0.0,
            effective_mass_kg=0.25 * ELECTRON_MASS_KG,
            dimension=dimension,
            degeneracy=2.0,
        )
        normalized = dos / dos[reference_index]
        axis.plot(
            excess,
            normalized,
            color=colors[dimension],
            lw=2,
            label=rf"{dimension}D: $g(E)\propto(E-E_0)^{{{dimension / 2 - 1:g}}}$",
        )
    axis.set_ylim(0.0, 5.0)
    axis.set_xlabel(r"energy above the band edge $E-E_0$ (eV)")
    axis.set_ylabel("DOS normalized at 0.25 eV")
    axis.set_title("Dimensionality changes the density of states")
    axis.legend(frameon=False)
    axis.grid(alpha=0.2)
    save(fig, "dos-dimensionality.pdf")


def numerical_convergence() -> None:
    """Expose the two numerical resolutions in a sampled-band DOS.

    Three panels independently report the number of reciprocal-space samples
    ``N_k`` and Gaussian width ``sigma``.  The sampled DOS is compared with the
    analytic cosine-band result, so mesh artifacts and excessive broadening are
    distinguishable.  State weights retain two spin-resolved states per cell
    in every panel.

    Returns
    -------
    None
        Writes ``dos-numerical-convergence.pdf`` into :data:`FIGURE_DIR`.
    """
    lattice_constant = 0.5e-9
    energy_grid = np.linspace(-2.4, 2.4, 1800)
    analytic_grid = np.linspace(-1.99, 1.99, 1500)
    analytic = tight_binding_dos_1d_analytic(
        analytic_grid,
        hopping_ev=1.0,
        degeneracy=2.0,
    )

    cases = [
        (64, 0.08, "coarse mesh"),
        (1024, 0.08, "converged mesh, broad kernel"),
        (8192, 0.025, "fine mesh, narrow kernel"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(9.4, 3.0), sharey=True, layout="constrained")
    for axis, (points, sigma, title) in zip(axes, cases, strict=True):
        wave_vector = uniform_bz_grid_1d(lattice_constant, points)
        band = tight_binding_band_1d(
            wave_vector,
            lattice_constant_m=lattice_constant,
            hopping_ev=1.0,
        )
        dos = sampled_band_density_of_states(
            band,
            energy_grid,
            sigma_ev=sigma,
            band_degeneracies=2.0,
        )
        axis.plot(energy_grid, dos, color="#123B66", lw=1.6, label="sampled")
        axis.plot(analytic_grid, analytic, color="#A44A17", lw=1, ls="--", label="analytic")
        axis.set_title(f"{title}\n$N_k={points}$, $\\sigma={sigma:.3f}$ eV")
        axis.set_xlabel(r"energy $E$ (eV)")
        axis.set_ylim(0.0, 2.6)
        axis.grid(alpha=0.18)
    axes[0].set_ylabel(r"DOS (states eV$^{-1}$ cell$^{-1}$)")
    axes[-1].legend(frameon=False, loc="upper center")
    save(fig, "dos-numerical-convergence.pdf")


def main() -> None:
    """Regenerate every Chapter 3 DOS figure.

    Returns
    -------
    None
        Writes all three DOS PDFs into :data:`FIGURE_DIR`.
    """
    band_to_dos()
    dimensionality_comparison()
    numerical_convergence()


if __name__ == "__main__":
    main()
