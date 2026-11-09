"""Local NumPy implementation of the LA1 two-dimensional Bloch model.

The module is intentionally colocated with the assessment notebooks.  It
constructs a plane-wave Hamiltonian for the checkerboard potential

    V(x, y) = V0 cos(2 pi x/a) cos(2 pi y/a)

in electron-volts and angstrom-based reciprocal coordinates.  It also provides
high-symmetry paths, endpoint-excluded Brillouin-zone meshes, Gaussian density
of states estimates, and stable Fermi--Dirac occupation.  No Project Koios or
PhysKit package is required.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.constants import Boltzmann, electron_mass, electron_volt, hbar

FloatArray = NDArray[np.float64]

HBAR2_OVER_2M0_EV_A2 = (hbar**2 / (2.0 * electron_mass)) / electron_volt / 1.0e-20
HBAR2_OVER_M0_EV_A2 = 2.0 * HBAR2_OVER_2M0_EV_A2
K_B_EV = Boltzmann / electron_volt


@dataclass(frozen=True, slots=True)
class CheckerboardBloch2D:
    """Plane-wave Bloch Hamiltonian for a square checkerboard potential.

    Parameters
    ----------
    lattice_constant_angstrom:
        Square-cell lattice constant in angstroms.
    potential_amplitude_ev:
        Coefficient ``V0`` in the checkerboard potential in electron-volts.
    reciprocal_mode_limit:
        Integer plane-wave cutoff ``-M <= m_x,m_y <= M``.
    """

    lattice_constant_angstrom: float = 5.43
    potential_amplitude_ev: float = 4.0
    reciprocal_mode_limit: int = 3

    def __post_init__(self) -> None:
        if self.lattice_constant_angstrom <= 0.0:
            raise ValueError("lattice_constant_angstrom must be positive")
        if not np.isfinite(self.potential_amplitude_ev):
            raise ValueError("potential_amplitude_ev must be finite")
        if type(self.reciprocal_mode_limit) is not int or self.reciprocal_mode_limit < 1:
            raise ValueError("reciprocal_mode_limit must be a positive integer")

    @property
    def zone_edge(self) -> float:
        """Return ``pi/a`` in inverse angstroms."""
        return float(np.pi / self.lattice_constant_angstrom)

    @property
    def zone_width(self) -> float:
        """Return the reciprocal width ``2 pi/a`` in inverse angstroms."""
        return 2.0 * self.zone_edge

    @property
    def direct_basis(self) -> FloatArray:
        """Return square direct-lattice basis vectors as matrix columns."""
        return self.lattice_constant_angstrom * np.eye(2)

    @property
    def reciprocal_basis(self) -> FloatArray:
        """Return reciprocal basis vectors as matrix columns."""
        return 2.0 * np.pi * np.linalg.inv(self.direct_basis).T

    @staticmethod
    def reciprocal_indices(mode_limit: int) -> NDArray[np.int64]:
        """Return integer reciprocal indices in deterministic lexicographic order."""
        if type(mode_limit) is not int or mode_limit < 1:
            raise ValueError("mode_limit must be a positive integer")
        return np.array(
            [
                (mx, my)
                for mx in range(-mode_limit, mode_limit + 1)
                for my in range(-mode_limit, mode_limit + 1)
            ],
            dtype=np.int64,
        )

    def fold_to_zone(self, wave_vector: ArrayLike) -> FloatArray:
        """Fold Cartesian wave vectors into ``[-pi/a, pi/a)``."""
        vector = np.asarray(wave_vector, dtype=float)
        return (vector + self.zone_edge) % self.zone_width - self.zone_edge

    def hamiltonian(
        self,
        wave_vector_inv_angstrom: ArrayLike,
        *,
        mode_limit: int | None = None,
    ) -> FloatArray:
        """Return the real symmetric plane-wave Hamiltonian in electron-volts."""
        k = np.asarray(wave_vector_inv_angstrom, dtype=float)
        if k.shape != (2,) or not np.all(np.isfinite(k)):
            raise ValueError("wave_vector_inv_angstrom must be a finite 2-vector")
        cutoff = self.reciprocal_mode_limit if mode_limit is None else mode_limit
        indices = self.reciprocal_indices(cutoff)
        reciprocal_vectors = indices @ self.reciprocal_basis.T
        kinetic = HBAR2_OVER_2M0_EV_A2 * np.sum(
            (reciprocal_vectors + k[None, :]) ** 2,
            axis=1,
        )
        matrix = np.diag(kinetic)
        differences = indices[:, None, :] - indices[None, :, :]
        coupled = (
            (np.abs(differences[:, :, 0]) == 1)
            & (np.abs(differences[:, :, 1]) == 1)
        )
        matrix[coupled] = self.potential_amplitude_ev / 4.0
        return matrix

    def solve_bands(
        self,
        wave_vectors_inv_angstrom: ArrayLike,
        *,
        band_count: int,
        mode_limit: int | None = None,
    ) -> FloatArray:
        """Diagonalize the Hamiltonian at each wave vector and return low bands."""
        vectors = np.atleast_2d(np.asarray(wave_vectors_inv_angstrom, dtype=float))
        if vectors.ndim != 2 or vectors.shape[1] != 2:
            raise ValueError("wave_vectors_inv_angstrom must have shape (N, 2)")
        cutoff = self.reciprocal_mode_limit if mode_limit is None else mode_limit
        basis_size = (2 * cutoff + 1) ** 2
        if type(band_count) is not int or not 1 <= band_count <= basis_size:
            raise ValueError("band_count must fit inside the represented basis")
        return np.array(
            [
                np.linalg.eigvalsh(
                    self.hamiltonian(vector, mode_limit=cutoff)
                )[:band_count]
                for vector in vectors
            ]
        )

    def high_symmetry_path(
        self,
        *,
        points_per_segment: int = 61,
    ) -> tuple[FloatArray, FloatArray, FloatArray, tuple[str, ...]]:
        """Return the square-lattice ``Gamma-X-M-Gamma`` path."""
        if type(points_per_segment) is not int or points_per_segment < 2:
            raise ValueError("points_per_segment must be at least two")
        z = self.zone_edge
        nodes = np.array(((0.0, 0.0), (z, 0.0), (z, z), (0.0, 0.0)))
        pieces = []
        for index, (start, stop) in enumerate(zip(nodes[:-1], nodes[1:], strict=True)):
            segment = np.linspace(start, stop, points_per_segment)
            if index:
                segment = segment[1:]
            pieces.append(segment)
        values = np.vstack(pieces)
        steps = np.linalg.norm(np.diff(values, axis=0), axis=1)
        distances = np.concatenate(([0.0], np.cumsum(steps)))
        node_distances = np.array(
            (
                0.0,
                np.linalg.norm(nodes[1] - nodes[0]),
                np.linalg.norm(nodes[1] - nodes[0])
                + np.linalg.norm(nodes[2] - nodes[1]),
                np.linalg.norm(nodes[1] - nodes[0])
                + np.linalg.norm(nodes[2] - nodes[1])
                + np.linalg.norm(nodes[3] - nodes[2]),
            )
        )
        return values, distances, node_distances, (r"$\Gamma$", r"$X$", r"$M$", r"$\Gamma$")

    def uniform_zone(
        self,
        axis_count: int,
        *,
        cell_centered: bool = False,
    ) -> tuple[FloatArray, FloatArray]:
        """Return one axis and all equal-weight endpoint-excluded zone vectors."""
        if type(axis_count) is not int or axis_count < 2:
            raise ValueError("axis_count must be at least two")
        if cell_centered:
            axis = (
                -self.zone_edge
                + (np.arange(axis_count) + 0.5) * self.zone_width / axis_count
            )
        else:
            axis = np.linspace(
                -self.zone_edge,
                self.zone_edge,
                axis_count,
                endpoint=False,
            )
        kx, ky = np.meshgrid(axis, axis, indexing="ij")
        return axis, np.column_stack((kx.ravel(), ky.ravel()))


def gaussian_dos_2d(
    band_energies_ev: ArrayLike,
    energy_grid_ev: ArrayLike,
    broadening_ev: float,
    sample_state_density_inv_angstrom2: float,
) -> FloatArray:
    """Return a normalized Gaussian DOS in states/(eV angstrom^2).

    ``sample_state_density_inv_angstrom2`` is the state weight assigned to one
    represented eigenvalue, including spin degeneracy.
    """
    energies = np.asarray(band_energies_ev, dtype=float).ravel()
    grid = np.asarray(energy_grid_ev, dtype=float)
    if energies.size == 0 or grid.ndim != 1 or grid.size < 2:
        raise ValueError("energies and a one-dimensional energy grid are required")
    if broadening_ev <= 0.0 or sample_state_density_inv_angstrom2 <= 0.0:
        raise ValueError("broadening and sample state density must be positive")
    result = np.zeros_like(grid)
    normalization = broadening_ev * np.sqrt(2.0 * np.pi)
    chunk_size = 512
    for start in range(0, energies.size, chunk_size):
        block = energies[start : start + chunk_size]
        reduced = (grid[:, None] - block[None, :]) / broadening_ev
        result += np.sum(np.exp(-0.5 * reduced**2), axis=1)
    return sample_state_density_inv_angstrom2 * result / normalization


def stable_fermi(
    energy_ev: ArrayLike,
    chemical_potential_ev: float,
    temperature_k: float,
) -> FloatArray:
    """Return overflow-safe Fermi--Dirac occupation."""
    if temperature_k <= 0.0:
        raise ValueError("temperature_k must be positive")
    argument = (
        np.asarray(energy_ev, dtype=float) - chemical_potential_ev
    ) / (K_B_EV * temperature_k)
    result = np.empty_like(argument)
    positive = argument >= 0.0
    exponential = np.exp(-argument[positive])
    result[positive] = exponential / (1.0 + exponential)
    exponential = np.exp(argument[~positive])
    result[~positive] = 1.0 / (1.0 + exponential)
    return result
