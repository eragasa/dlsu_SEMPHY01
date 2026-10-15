"""Reproducible figures for the particle-box-to-band sequence.

This module is the reusable implementation behind the student notebooks for
Chapters 1 and 2.  Its figures connect three descriptions of one-electron
states: standing waves in a hard-wall box, translation eigenstates in a
periodic crystal, and molecular orbitals produced by a localized tight-binding
basis.  Energies in the schematic models are dimensionless or normalized by
the hopping amplitude, as stated on each axis.

The plotting functions accept no arguments because each published figure is a
fixed scientific artifact.  Every numerical value appearing here is therefore
part of the figure's reproducibility record rather than an implicit notebook
state.  Each function writes one PDF into :data:`FIGURE_DIR` and closes its
Matplotlib figure.  The PDF metadata omits creation and modification times so
that repeated executions produce byte-identical files.

Students should first work through the direct calculation in the corresponding
``figures/0X_chapter-name/fig_0X_0X.ipynb`` notebook.  The final notebook cell
imports the appropriate function from this module to demonstrate how the same
calculation is packaged for reliable reuse.
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
        "figure.dpi": 120,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
        "pdf.fonttype": 42,
    }
)

METADATA = {
    "Creator": "SEMPHY01 deterministic band-foundation figure generator",
    "Producer": "Matplotlib",
    "CreationDate": None,
    "ModDate": None,
}


def save(figure: plt.Figure, name: str) -> None:
    """Write and close one deterministic publication figure.

    Parameters
    ----------
    figure:
        Fully constructed Matplotlib figure.  This function takes ownership of
        the figure and closes it after writing, preventing hidden state from
        leaking into the next notebook or test.
    name:
        PDF filename relative to :data:`FIGURE_DIR`.  Callers pass a basename,
        keeping every generated PDF directly under ``figures/`` where LaTeX
        can include it with a stable path.

    Returns
    -------
    None
        The observable result is the PDF written to disk.

    Notes
    -----
    Fixed metadata removes timestamp variation.  Deterministic numerical
    arrays and fixed style parameters then make the output suitable for hash
    comparison in reproducibility checks.
    """
    figure.savefig(FIGURE_DIR / name, metadata=METADATA)
    plt.close(figure)


def piab_spectrum() -> None:
    """Plot normalized hard-wall eigenfunctions and their discrete energies.

    The calculation uses the dimensionless coordinate ``x/L``.  The left panel
    evaluates ``sqrt(2) sin(n pi x/L)`` for ``n = 1, ..., 4`` and offsets the
    curves vertically so their node structure remains visible.  The right
    panel plots ``E_n/E_1 = n^2``, displaying the quadratic separation of the
    box spectrum.

    Returns
    -------
    None
        Writes ``foundation-piab-spectrum.pdf`` into :data:`FIGURE_DIR`.
    """
    x = np.linspace(0.0, 1.0, 1000)
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8), layout="constrained")

    for n in range(1, 5):
        wave = np.sqrt(2.0) * np.sin(n * np.pi * x)
        axes[0].plot(x, wave + 3.0 * (n - 1), lw=1.7, label=rf"$n={n}$")
        axes[0].axhline(3.0 * (n - 1), color="0.75", lw=0.6)
    axes[0].set_xlabel(r"position $x/L$")
    axes[0].set_ylabel("wavefunction, vertically offset")
    axes[0].set_title("standing-wave eigenfunctions")
    axes[0].set_yticks([])
    axes[0].legend(frameon=False, ncol=2)

    for n in range(1, 5):
        axes[1].hlines(n**2, 0.15, 0.85, color="#123B66", lw=2)
        axes[1].text(0.89, n**2, rf"$n={n}$", va="center")
    axes[1].set_xlim(0.0, 1.1)
    axes[1].set_ylim(0.0, 17.5)
    axes[1].set_xticks([])
    axes[1].set_ylabel(r"energy $E_n/E_1=n^2$")
    axes[1].set_title("discrete energy spectrum")
    save(fig, "foundation-piab-spectrum.pdf")


def boundary_conditions() -> None:
    """Contrast hard-wall standing waves with periodic plane waves.

    The hard-wall example is the third sine eigenfunction, whose amplitude
    vanishes at both ends of the interval.  The periodic example is the complex
    plane wave ``exp(i 4 pi x/L)``; separate real and imaginary curves expose
    its phase while endpoint markers verify periodicity.

    Returns
    -------
    None
        Writes ``foundation-boundary-conditions.pdf`` into
        :data:`FIGURE_DIR`.
    """
    x = np.linspace(0.0, 1.0, 1000)
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.2), layout="constrained")

    standing = np.sqrt(2.0) * np.sin(3.0 * np.pi * x)
    axes[0].plot(x, standing, color="#123B66", lw=2)
    axes[0].scatter([0.0, 1.0], [0.0, 0.0], color="#A44A17", zorder=3)
    axes[0].axhline(0.0, color="0.6", lw=0.7)
    axes[0].set_title(r"hard walls: $k_3=3\pi/L$")
    axes[0].set_xlabel(r"$x/L$")
    axes[0].set_ylabel(r"$\phi(x)$")

    real = np.cos(4.0 * np.pi * x)
    imag = np.sin(4.0 * np.pi * x)
    axes[1].plot(x, real, color="#123B66", lw=2, label="real part")
    axes[1].plot(x, imag, color="#A44A17", lw=1.6, ls="--", label="imaginary part")
    axes[1].scatter([0.0, 1.0], [1.0, 1.0], color="#2E7D32", zorder=3)
    axes[1].axhline(0.0, color="0.6", lw=0.7)
    axes[1].set_title(r"periodic box: $k_2=4\pi/L$")
    axes[1].set_xlabel(r"$x/L$")
    axes[1].set_ylabel(r"$e^{ikx}$")
    axes[1].legend(frameon=False)
    save(fig, "foundation-boundary-conditions.pdf")


def bloch_wave() -> None:
    """Display the factors in the one-dimensional Bloch form.

    A positive lattice-periodic function ``u_k(x)`` with period ``a`` is
    multiplied by the real part of ``exp(i k x)``.  Three aligned panels show
    the periodic factor, the long-wavelength phase, and their product
    ``Re[psi_k(x)]`` across six cells.  Vertical guides identify the lattice
    translations.

    Returns
    -------
    None
        Writes ``foundation-bloch-wave.pdf`` into :data:`FIGURE_DIR`.
    """
    cells = 6
    x = np.linspace(0.0, float(cells), 1800)
    periodic = 1.0 + 0.38 * np.cos(2.0 * np.pi * x)
    phase = np.cos(0.55 * np.pi * x)
    bloch_real = periodic * phase

    fig, axes = plt.subplots(3, 1, figsize=(8.2, 4.5), sharex=True, layout="constrained")
    axes[0].plot(x, periodic, color="#2E7D32", lw=1.8)
    axes[0].set_ylabel(r"$u_k(x)$")
    axes[0].set_title("lattice-periodic factor")
    axes[1].plot(x, phase, color="#A44A17", lw=1.8)
    axes[1].set_ylabel(r"$\Re\,e^{ikx}$")
    axes[1].set_title("plane-wave phase")
    axes[2].plot(x, bloch_real, color="#123B66", lw=2)
    axes[2].set_ylabel(r"$\Re\,\psi_k(x)$")
    axes[2].set_xlabel(r"position $x/a$")
    axes[2].set_title("Bloch state")
    for axis in axes:
        axis.axhline(0.0, color="0.7", lw=0.6)
        for boundary in range(cells + 1):
            axis.axvline(boundary, color="0.85", lw=0.5)
    save(fig, "foundation-bloch-wave.pdf")


def two_site_splitting() -> None:
    """Draw the bonding-antibonding splitting of two equivalent orbitals.

    The two-site Hamiltonian has on-site energy ``epsilon`` and hopping
    matrix element ``-t``.  Its symmetric and antisymmetric eigenstates have
    energies ``epsilon - t`` and ``epsilon + t``; their separation is ``2t``.
    The energy diagram makes this eigenvalue calculation visible.

    Returns
    -------
    None
        Writes ``foundation-two-site-splitting.pdf`` into
        :data:`FIGURE_DIR`.
    """
    epsilon = 0.0
    hopping = 1.0
    fig, axis = plt.subplots(figsize=(6.8, 3.5), layout="constrained")
    axis.hlines(epsilon, 0.08, 0.30, color="#666666", lw=2)
    axis.hlines(epsilon, 0.70, 0.92, color="#666666", lw=2)
    axis.hlines(epsilon - hopping, 0.38, 0.62, color="#123B66", lw=2.5)
    axis.hlines(epsilon + hopping, 0.38, 0.62, color="#A44A17", lw=2.5)
    axis.plot([0.30, 0.38], [epsilon, epsilon - hopping], color="0.6", lw=0.8)
    axis.plot([0.70, 0.62], [epsilon, epsilon - hopping], color="0.6", lw=0.8)
    axis.plot([0.30, 0.38], [epsilon, epsilon + hopping], color="0.6", lw=0.8)
    axis.plot([0.70, 0.62], [epsilon, epsilon + hopping], color="0.6", lw=0.8)
    axis.text(0.19, epsilon + 0.10, r"$|A\rangle$", ha="center")
    axis.text(0.81, epsilon + 0.10, r"$|B\rangle$", ha="center")
    axis.text(0.64, epsilon - hopping, r"$\epsilon-t$ bonding", va="center")
    axis.text(0.64, epsilon + hopping, r"$\epsilon+t$ antibonding", va="center")
    axis.annotate(
        "",
        xy=(0.50, epsilon + hopping),
        xytext=(0.50, epsilon - hopping),
        arrowprops={"arrowstyle": "<->", "color": "#2E7D32"},
    )
    axis.text(0.47, epsilon, r"splitting $2t$", ha="right", va="center")
    axis.set_xlim(0.0, 1.0)
    axis.set_ylim(-1.5, 1.5)
    axis.set_xticks([])
    axis.set_ylabel(r"energy relative to $\epsilon$")
    axis.set_title("two-site molecular-orbital splitting")
    axis.spines[["top", "right", "bottom"]].set_visible(False)
    save(fig, "foundation-two-site-splitting.pdf")


def levels_to_band() -> None:
    """Show finite-chain eigenvalues becoming dense inside a band.

    For an open chain of ``N`` sites with zero on-site energy and unit hopping,
    the plotted levels are ``E_m/t = -2 cos[m pi/(N+1)]``.  Columns for
    ``N = 2, 5, 12, 50`` preserve the exact count of one orbital per site and
    approach the limiting interval from ``-2t`` to ``+2t``.

    Returns
    -------
    None
        Writes ``foundation-levels-to-band.pdf`` into :data:`FIGURE_DIR`.
    """
    site_counts = [2, 5, 12, 50]
    fig, axis = plt.subplots(figsize=(7.4, 3.8), layout="constrained")
    for column, count in enumerate(site_counts):
        mode = np.arange(1, count + 1)
        energies = -2.0 * np.cos(mode * np.pi / (count + 1))
        width = 0.24 if count < 20 else 0.32
        axis.hlines(
            energies,
            column - width,
            column + width,
            color="#123B66",
            lw=1.2 if count < 20 else 0.55,
        )
    axis.axhline(-2.0, color="#A44A17", ls="--", lw=0.9)
    axis.axhline(2.0, color="#A44A17", ls="--", lw=0.9)
    axis.set_xticks(range(len(site_counts)), [rf"$N={n}$" for n in site_counts])
    axis.set_ylabel(r"$(E-\epsilon)/t$")
    axis.set_title("finite molecular levels approaching a band")
    axis.set_ylim(-2.3, 2.3)
    axis.grid(alpha=0.15, axis="y")
    save(fig, "foundation-levels-to-band.pdf")


def dimerized_bands() -> None:
    """Compare gapless and gapped two-orbital tight-binding bands.

    The Bloch eigenvalues are ``E_±(k) = ±|t1 + t2 exp(-ika)|``.  Equal
    intracell and intercell couplings close the gap at the zone boundary;
    unequal couplings open the analytic gap ``2|t1-t2|``.  The two panels use
    the same energy scale so this change is quantitative rather than merely
    schematic.

    Returns
    -------
    None
        Writes ``foundation-dimerized-bands.pdf`` into :data:`FIGURE_DIR`.
    """
    reduced_k = np.linspace(-1.0, 1.0, 1200)
    ka = np.pi * reduced_k
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.4), sharey=True, layout="constrained")
    cases = [(1.0, 1.0, "equal couplings"), (1.3, 0.7, "alternating couplings")]
    for axis, (t1, t2, title) in zip(axes, cases, strict=True):
        magnitude = np.sqrt(t1**2 + t2**2 + 2.0 * t1 * t2 * np.cos(ka))
        axis.plot(reduced_k, -magnitude, color="#123B66", lw=2, label=r"$E_-(k)$")
        axis.plot(reduced_k, magnitude, color="#A44A17", lw=2, label=r"$E_+(k)$")
        axis.axvline(-1.0, color="0.65", ls=":", lw=0.7)
        axis.axvline(1.0, color="0.65", ls=":", lw=0.7)
        axis.set_xlabel(r"$ka/\pi$")
        axis.set_title(title)
        axis.grid(alpha=0.15)
    axes[0].set_ylabel(r"energy $(E-\epsilon)$")
    gap_half = abs(cases[1][0] - cases[1][1])
    axes[1].annotate(
        "",
        xy=(0.92, gap_half),
        xytext=(0.92, -gap_half),
        arrowprops={"arrowstyle": "<->", "color": "#2E7D32"},
    )
    axes[1].text(0.48, 0.0, r"gap $2|t_1-t_2|$", ha="center", va="center")
    axes[0].legend(frameon=False, loc="upper center")
    save(fig, "foundation-dimerized-bands.pdf")


def main() -> None:
    """Regenerate every fixed figure provided by this module.

    This convenience entry point supports manual verification.  The manuscript
    workflow ordinarily invokes each recipe from its one-figure notebook so
    the explanatory calculation and the published artifact remain adjacent.

    Returns
    -------
    None
        Writes all six foundation PDFs into :data:`FIGURE_DIR`.
    """
    piab_spectrum()
    boundary_conditions()
    bloch_wave()
    two_site_splitting()
    levels_to_band()
    dimerized_bands()


if __name__ == "__main__":
    main()
