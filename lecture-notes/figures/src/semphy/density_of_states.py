"""Reusable numerical relations for bands and densities of states.

This module contains the scientific calculations shared by the Chapter 3
notebooks and publication figures.  The central distinction is between a
finite list of sampled eigenvalues and a continuum density of states.  Sampled
calculations attach an explicit state weight to every eigenvalue.  Continuum
calculations derive the density from reciprocal-space volume and therefore
carry a spatial-density unit.

Unit conventions
----------------
Energies are expressed in electronvolts, wave vectors in inverse metres,
lengths in metres, and effective masses in kilograms.  Functions that return a
continuum parabolic DOS use states per electronvolt per metre raised to the
requested spatial dimension.  Functions for a sampled band return states per
electronvolt per modeled primitive cell because their caller-supplied weights
encode the number of states represented by each sample.

Normalization conventions
-------------------------
A single spin-degenerate band contains two states per primitive cell.  For
``N_k`` equal-weight samples, each energy therefore receives weight
``2/N_k``.  Multiband calculations sum the declared degeneracy of every band.
The module never inserts an implicit valley or orbital degeneracy; callers must
state each such factor explicitly.

Numerical scope
---------------
The routines implement educational, inspectable calculations.  Gaussian
broadening is a representation of a discrete sample, not a physical lifetime.
A finite reciprocal-space grid controls numerical resolution, while an actual
spectral width would require an additional physical model.  Input validation
fails early on nonfinite values, inconsistent shapes, invalid dimensions, and
ambiguous normalizations.
"""

from __future__ import annotations

from collections.abc import Sequence
from math import gamma

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]

ELEMENTARY_CHARGE_J_PER_EV = 1.602176634e-19
HBAR_J_S = 1.054571817e-34
ELECTRON_MASS_KG = 9.1093837139e-31


def _float_array(values: ArrayLike, *, name: str) -> FloatArray:
    """Convert numerical input to a finite double-precision NumPy array.

    Parameters
    ----------
    values:
        Scalar or array-like numerical input accepted by ``numpy.asarray``.
    name:
        Physical variable name used in a diagnostic if a NaN or infinity is
        found.

    Returns
    -------
    numpy.ndarray
        A ``float64`` view or copy with the same shape as ``values``.

    Raises
    ------
    ValueError
        If any element is nonfinite.  Silent propagation of NaNs would make a
        DOS normalization check physically meaningless, so validation occurs
        at the computational boundary.
    """
    array = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def uniform_bz_grid_1d(lattice_constant_m: float, points: int) -> FloatArray:
    """Return an equal-weight midpoint grid for the first one-dimensional BZ.

    Parameters
    ----------
    lattice_constant_m:
        Positive primitive-cell length ``a`` in metres.
    points:
        Number ``N_k`` of reciprocal-space samples.  At least two points are
        required.

    Returns
    -------
    numpy.ndarray
        Wave vectors in inverse metres.  The points lie at cell midpoints over
        ``[-pi/a, pi/a)`` and each represents the Brillouin-zone fraction
        ``1/N_k``.

    Raises
    ------
    ValueError
        If ``a`` is nonpositive or fewer than two points are requested.

    Notes
    -----
    Midpoint sampling avoids duplicating the equivalent boundaries ``-pi/a``
    and ``+pi/a``.  It is particularly convenient for state-count sums because
    all samples have identical weights.
    """
    if lattice_constant_m <= 0.0:
        raise ValueError("lattice_constant_m must be positive")
    if points < 2:
        raise ValueError("points must be at least two")
    spacing = 2.0 * np.pi / (lattice_constant_m * points)
    return -np.pi / lattice_constant_m + (np.arange(points) + 0.5) * spacing


def tight_binding_band_1d(
    wave_vector_per_m: ArrayLike,
    *,
    lattice_constant_m: float,
    onsite_ev: float = 0.0,
    hopping_ev: float = 1.0,
) -> FloatArray:
    """Evaluate the nearest-neighbor one-orbital Bloch band.

    Parameters
    ----------
    wave_vector_per_m:
        Wave vector or array of wave vectors ``k`` in inverse metres.
    lattice_constant_m:
        Positive lattice constant ``a`` in metres.
    onsite_ev:
        Localized-orbital on-site energy ``epsilon`` in electronvolts.
    hopping_ev:
        Positive hopping magnitude ``t`` in electronvolts.  The Hamiltonian
        convention uses matrix elements ``-t`` between neighbors.

    Returns
    -------
    numpy.ndarray
        Energies ``epsilon - 2t cos(ka)`` in electronvolts with the same shape
        as the input wave-vector array.

    Raises
    ------
    ValueError
        If the lattice constant or hopping is nonpositive, or if the wave
        vectors contain NaN or infinity.

    Notes
    -----
    This is one solvable model for a Bloch band.  It is not the definition of a
    band and carries no occupation information.
    """
    if lattice_constant_m <= 0.0:
        raise ValueError("lattice_constant_m must be positive")
    if hopping_ev <= 0.0:
        raise ValueError("hopping_ev must be positive")
    wave_vector = _float_array(wave_vector_per_m, name="wave_vector_per_m")
    return onsite_ev - 2.0 * hopping_ev * np.cos(
        wave_vector * lattice_constant_m
    )


def tight_binding_dos_1d_analytic(
    energy_ev: ArrayLike,
    *,
    onsite_ev: float = 0.0,
    hopping_ev: float = 1.0,
    degeneracy: float = 2.0,
) -> FloatArray:
    """Return the analytic DOS of a nearest-neighbor cosine band.

    Parameters
    ----------
    energy_ev:
        Energies at which to evaluate the density of states, in electronvolts.
    onsite_ev:
        Band-center energy ``epsilon`` in electronvolts.
    hopping_ev:
        Positive hopping magnitude ``t`` in electronvolts.  The allowed band
        extends from ``epsilon-2t`` to ``epsilon+2t``.
    degeneracy:
        Number of states represented by each spatial Bloch state.  Use two for
        a spin-degenerate band and one for a spin-resolved band.

    Returns
    -------
    numpy.ndarray
        DOS in states per electronvolt per primitive cell.  Values outside the
        band are zero, values inside follow the inverse square-root law, and
        values exactly at an ideal band edge are infinite.

    Raises
    ------
    ValueError
        If hopping or degeneracy is nonpositive, or energies are nonfinite.

    Notes
    -----
    Integrating across the complete band gives ``degeneracy``.  The divergent
    edges are one-dimensional van Hove singularities produced by a vanishing
    band slope; a finite numerical broadening rounds them.
    """
    if hopping_ev <= 0.0:
        raise ValueError("hopping_ev must be positive")
    if degeneracy <= 0.0:
        raise ValueError("degeneracy must be positive")
    energy = _float_array(energy_ev, name="energy_ev")
    radicand = (2.0 * hopping_ev) ** 2 - (energy - onsite_ev) ** 2
    result = np.zeros_like(energy)
    inside = radicand > 0.0
    result[inside] = degeneracy / (np.pi * np.sqrt(radicand[inside]))
    at_edge = np.isclose(radicand, 0.0, rtol=0.0, atol=1.0e-14)
    result[at_edge] = np.inf
    return result


def gaussian_density_of_states(
    energies_ev: ArrayLike,
    energy_grid_ev: ArrayLike,
    *,
    sigma_ev: float,
    state_weights: ArrayLike | None = None,
    chunk_size: int = 4096,
) -> FloatArray:
    """Represent weighted discrete eigenvalues by normalized Gaussian kernels.

    Parameters
    ----------
    energies_ev:
        One-dimensional collection of discrete energies in electronvolts.
    energy_grid_ev:
        Strictly increasing output grid in electronvolts.
    sigma_ev:
        Positive Gaussian standard deviation in electronvolts.  This parameter
        controls numerical smoothing; it has no automatic interpretation as a
        physical lifetime.
    state_weights:
        Optional nonnegative number of physical states represented by each
        energy.  If omitted, every eigenvalue has unit weight.
    chunk_size:
        Maximum number of kernels evaluated simultaneously.  Chunking bounds
        temporary memory without changing the mathematical result.

    Returns
    -------
    numpy.ndarray
        Gaussian-broadened DOS evaluated on ``energy_grid_ev``, in the state
        unit supplied by ``state_weights`` per electronvolt.

    Raises
    ------
    ValueError
        If an input is empty or nonfinite, the grid is not increasing, sigma or
        chunk size is nonpositive, or weights have the wrong length or sign.

    Notes
    -----
    Every kernel integrates to its state weight on an infinite energy domain.
    A finite output interval must extend several ``sigma_ev`` beyond the
    extreme eigenvalues if numerical integration is expected to recover the
    complete state count.
    """
    energies = _float_array(energies_ev, name="energies_ev").reshape(-1)
    grid = _float_array(energy_grid_ev, name="energy_grid_ev").reshape(-1)
    if energies.size == 0 or grid.size < 2:
        raise ValueError("energies must be nonempty and grid must have two points")
    if not np.all(np.diff(grid) > 0.0):
        raise ValueError("energy_grid_ev must be strictly increasing")
    if sigma_ev <= 0.0:
        raise ValueError("sigma_ev must be positive")
    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")

    if state_weights is None:
        weights = np.ones(energies.size)
    else:
        weights = _float_array(state_weights, name="state_weights").reshape(-1)
        if weights.size != energies.size:
            raise ValueError("state_weights must match energies_ev")
        if np.any(weights < 0.0):
            raise ValueError("state_weights must be nonnegative")

    result = np.zeros_like(grid)
    normalization = sigma_ev * np.sqrt(2.0 * np.pi)
    for start in range(0, energies.size, chunk_size):
        stop = min(start + chunk_size, energies.size)
        reduced = (
            grid[:, np.newaxis] - energies[np.newaxis, start:stop]
        ) / sigma_ev
        kernels = np.exp(-0.5 * reduced**2) / normalization
        result += kernels @ weights[start:stop]
    return result


def histogram_density_of_states(
    energies_ev: ArrayLike,
    bin_edges_ev: ArrayLike,
    *,
    state_weights: ArrayLike | None = None,
) -> tuple[FloatArray, FloatArray]:
    """Convert weighted eigenvalues into a piecewise-constant histogram DOS.

    Parameters
    ----------
    energies_ev:
        Discrete energy samples in electronvolts.
    bin_edges_ev:
        Strictly increasing histogram boundaries in electronvolts.
    state_weights:
        Optional weight carried by each energy.  Unit weights count individual
        discrete states; Brillouin-zone weights produce states per cell.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        Bin centers in electronvolts and DOS values equal to weighted counts
        divided by each bin width.

    Raises
    ------
    ValueError
        If edges are invalid, values are nonfinite, or the weight array has the
        wrong length.

    Notes
    -----
    Multiplying the returned DOS by the corresponding bin widths and summing
    recovers the weight of all energies that fall inside the supplied range.
    """
    energies = _float_array(energies_ev, name="energies_ev").reshape(-1)
    edges = _float_array(bin_edges_ev, name="bin_edges_ev").reshape(-1)
    if edges.size < 2 or not np.all(np.diff(edges) > 0.0):
        raise ValueError("bin_edges_ev must be strictly increasing")
    if state_weights is None:
        weights = np.ones(energies.size)
    else:
        weights = _float_array(state_weights, name="state_weights").reshape(-1)
        if weights.size != energies.size:
            raise ValueError("state_weights must match energies_ev")
    counts, _ = np.histogram(energies, bins=edges, weights=weights)
    widths = np.diff(edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, counts / widths


def sampled_band_density_of_states(
    band_energies_ev: ArrayLike,
    energy_grid_ev: ArrayLike,
    *,
    sigma_ev: float,
    k_weights: ArrayLike | None = None,
    band_degeneracies: float | Sequence[float] = 2.0,
) -> FloatArray:
    """Calculate a normalized Gaussian DOS from one or more sampled bands.

    Parameters
    ----------
    band_energies_ev:
        Energy array with shape ``(bands, k_points)`` or ``(k_points,)`` in
        electronvolts.  A one-dimensional input represents one band.
    energy_grid_ev:
        Strictly increasing output grid in electronvolts.
    sigma_ev:
        Positive Gaussian broadening width in electronvolts.
    k_weights:
        Optional nonnegative integration weights for reciprocal-space points.
        They must sum to one.  Equal weights are supplied by default.
    band_degeneracies:
        Positive scalar applied to every band or one positive value per band.
        This is where spin, valley, or explicitly repeated band copies enter.

    Returns
    -------
    numpy.ndarray
        DOS on the requested energy grid in states per electronvolt per modeled
        cell.

    Raises
    ------
    ValueError
        If shapes, weights, degeneracies, or numerical values violate the
        stated contracts.

    Notes
    -----
    The integrated DOS equals the sum of ``band_degeneracies`` when the energy
    grid captures the complete broadened spectrum.  This state-count sum rule
    is the primary normalization test for a sampled band.
    """
    bands = _float_array(band_energies_ev, name="band_energies_ev")
    if bands.ndim == 1:
        bands = bands[np.newaxis, :]
    if bands.ndim != 2 or bands.shape[1] < 2:
        raise ValueError("band_energies_ev must have shape (bands, k_points)")
    band_count, k_count = bands.shape

    if k_weights is None:
        weights_k = np.full(k_count, 1.0 / k_count)
    else:
        weights_k = _float_array(k_weights, name="k_weights").reshape(-1)
        if weights_k.size != k_count or np.any(weights_k < 0.0):
            raise ValueError("k_weights must be nonnegative and match k_points")
        if not np.isclose(np.sum(weights_k), 1.0, rtol=1.0e-12, atol=1.0e-14):
            raise ValueError("k_weights must sum to one")

    degeneracies = np.asarray(band_degeneracies, dtype=float)
    if degeneracies.ndim == 0:
        degeneracies = np.full(band_count, float(degeneracies))
    if degeneracies.shape != (band_count,) or np.any(degeneracies <= 0.0):
        raise ValueError("band_degeneracies must be positive and match bands")

    state_weights = (degeneracies[:, np.newaxis] * weights_k).reshape(-1)
    return gaussian_density_of_states(
        bands.reshape(-1),
        energy_grid_ev,
        sigma_ev=sigma_ev,
        state_weights=state_weights,
    )


def parabolic_band_energy(
    wave_vectors_per_m: ArrayLike,
    *,
    edge_ev: float,
    principal_masses_kg: float | Sequence[float],
) -> FloatArray:
    """Evaluate an anisotropic parabolic band in its principal-axis basis.

    Parameters
    ----------
    wave_vectors_per_m:
        Wave-vector components in inverse metres.  The final axis enumerates
        Cartesian principal directions; earlier axes enumerate sample points.
    edge_ev:
        Band-edge energy in electronvolts.
    principal_masses_kg:
        One positive isotropic mass or one positive effective mass per final-
        axis component, in kilograms.

    Returns
    -------
    numpy.ndarray
        Energies ``edge + hbar^2 sum(k_i^2/m_i)/2`` in electronvolts, with the
        component axis removed.

    Raises
    ------
    ValueError
        If input values are nonfinite, masses are nonpositive, or the mass
        vector does not match the number of components.

    Notes
    -----
    The expression describes a local quadratic expansion of a band extremum.
    It does not assert that the complete band remains parabolic far from that
    extremum.
    """
    wave_vectors = _float_array(
        wave_vectors_per_m, name="wave_vectors_per_m"
    )
    if wave_vectors.ndim == 0:
        wave_vectors = wave_vectors.reshape(1, 1)
    elif wave_vectors.ndim == 1:
        wave_vectors = wave_vectors[:, np.newaxis]
    masses = np.asarray(principal_masses_kg, dtype=float)
    if masses.ndim == 0:
        masses = np.full(wave_vectors.shape[-1], float(masses))
    if masses.shape != (wave_vectors.shape[-1],) or np.any(masses <= 0.0):
        raise ValueError("principal masses must be positive and match dimension")
    kinetic_j = 0.5 * HBAR_J_S**2 * np.sum(
        wave_vectors**2 / masses, axis=-1
    )
    return edge_ev + kinetic_j / ELEMENTARY_CHARGE_J_PER_EV


def density_of_states_mass(principal_masses_kg: Sequence[float]) -> float:
    """Return the geometric-mean mass controlling an ellipsoidal DOS.

    Parameters
    ----------
    principal_masses_kg:
        Nonempty sequence of positive principal effective masses in kilograms.

    Returns
    -------
    float
        ``(m1*m2*...*md)^(1/d)`` in kilograms.

    Raises
    ------
    ValueError
        If the input is empty, multidimensional, or contains a nonpositive
        mass.

    Notes
    -----
    In three dimensions this definition satisfies
    ``m_dos^(3/2) = (m1*m2*m3)^(1/2)``, the factor appearing in the parabolic
    state count of an ellipsoidal constant-energy surface.
    """
    masses = np.asarray(principal_masses_kg, dtype=float)
    if masses.ndim != 1 or masses.size == 0 or np.any(masses <= 0.0):
        raise ValueError("principal_masses_kg must be a nonempty positive vector")
    return float(np.prod(masses) ** (1.0 / masses.size))


def parabolic_density_of_states(
    energy_ev: ArrayLike,
    *,
    edge_ev: float,
    effective_mass_kg: float,
    dimension: int,
    degeneracy: float = 2.0,
) -> FloatArray:
    """Return the continuum DOS of an isotropic parabolic band edge.

    Parameters
    ----------
    energy_ev:
        Evaluation energies in electronvolts.
    edge_ev:
        Band-edge energy in electronvolts.
    effective_mass_kg:
        Positive scalar effective mass in kilograms.  For an anisotropic
        ellipsoid, pass the DOS mass appropriate to the selected dimension.
    dimension:
        Spatial dimension ``1``, ``2``, or ``3``.
    degeneracy:
        Explicit spin, valley, or equivalent-band multiplicity.

    Returns
    -------
    numpy.ndarray
        DOS in states per electronvolt per metre raised to ``dimension``.
        Values below the edge are zero; the ideal 1D value at the edge is
        infinite.

    Raises
    ------
    ValueError
        If dimension is unsupported, mass or degeneracy is nonpositive, or an
        energy is nonfinite.

    Notes
    -----
    The derivation counts the reciprocal-space volume of a ``d``-dimensional
    sphere and differentiates with respect to energy.  It yields the universal
    edge law ``g(E) proportional to (E-edge)^(d/2-1)``.
    """
    if dimension not in (1, 2, 3):
        raise ValueError("dimension must be 1, 2, or 3")
    if effective_mass_kg <= 0.0 or degeneracy <= 0.0:
        raise ValueError("effective_mass_kg and degeneracy must be positive")
    energy = _float_array(energy_ev, name="energy_ev")
    excess_j = (energy - edge_ev) * ELEMENTARY_CHARGE_J_PER_EV
    result = np.zeros_like(energy)
    positive = excess_j > 0.0

    d = float(dimension)
    unit_ball_volume = np.pi ** (0.5 * d) / gamma(0.5 * d + 1.0)
    coefficient_per_joule = (
        degeneracy
        / (2.0 * np.pi) ** d
        * unit_ball_volume
        * (0.5 * d)
        * (2.0 * effective_mass_kg / HBAR_J_S**2) ** (0.5 * d)
    )
    result[positive] = (
        coefficient_per_joule
        * excess_j[positive] ** (0.5 * d - 1.0)
        * ELEMENTARY_CHARGE_J_PER_EV
    )
    if dimension == 1:
        result[np.isclose(excess_j, 0.0, rtol=0.0, atol=0.0)] = np.inf
    return result


def integrated_states(energy_ev: ArrayLike, dos_per_ev: ArrayLike) -> float:
    """Integrate a tabulated density of states over energy.

    Parameters
    ----------
    energy_ev:
        Strictly increasing energy grid in electronvolts.
    dos_per_ev:
        Nonnegative DOS values on the same grid.

    Returns
    -------
    float
        Trapezoidal estimate of the represented state count.  Its remaining
        spatial unit is inherited from the DOS.

    Raises
    ------
    ValueError
        If arrays have different lengths, contain fewer than two values, are
        nonfinite, use a nonmonotonic energy grid, or contain a negative DOS.

    Notes
    -----
    This function is intentionally small because normalization is a scientific
    invariant: a complete isolated band must integrate to its declared
    degeneracy.  Callers should refine both the energy grid and spectral range
    until this value is stable.
    """
    energy = _float_array(energy_ev, name="energy_ev").reshape(-1)
    dos = _float_array(dos_per_ev, name="dos_per_ev").reshape(-1)
    if energy.size != dos.size or energy.size < 2:
        raise ValueError("energy_ev and dos_per_ev must have matching lengths")
    if not np.all(np.diff(energy) > 0.0):
        raise ValueError("energy_ev must be strictly increasing")
    if np.any(dos < 0.0):
        raise ValueError("dos_per_ev must be nonnegative")
    return float(np.trapezoid(dos, energy))
