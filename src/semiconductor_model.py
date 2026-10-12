"""Reusable equilibrium-semiconductor model scaffold.

The public data structures are usable now; physics and solver functions marked
with ``NotImplementedError`` are explicit Part II/IV implementation tasks.
Transport quantities do not belong in this module.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True, slots=True)
class BandValley:
    """One set of equivalent parabolic 3D band-edge valleys."""

    label: str
    edge_energy_ev: float
    principal_masses_m0: tuple[float, float, float]
    valley_degeneracy: int = 1

    def __post_init__(self) -> None:
        if len(self.principal_masses_m0) != 3:
            raise ValueError("principal_masses_m0 must contain three masses")
        if any(mass <= 0.0 for mass in self.principal_masses_m0):
            raise ValueError("principal masses must be positive")
        if self.valley_degeneracy < 1:
            raise ValueError("valley_degeneracy must be positive")

    @property
    def dos_mass_m0(self) -> float:
        """Return the per-valley 3D DOS mass in free-electron-mass units."""
        return float(np.prod(self.principal_masses_m0) ** (1.0 / 3.0))


@dataclass(frozen=True, slots=True)
class BandEdgeModel:
    """Bulk band edges and all modeled conduction/valence components."""

    material: str
    phase_or_polytype: str
    valence_edge_ev: float
    conduction_edge_ev: float
    conduction_valleys: tuple[BandValley, ...]
    valence_bands: tuple[BandValley, ...]
    provenance: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.conduction_edge_ev <= self.valence_edge_ev:
            raise ValueError("conduction edge must lie above valence edge")
        if not self.conduction_valleys or not self.valence_bands:
            raise ValueError("at least one conduction and valence component is required")

    @property
    def gap_ev(self) -> float:
        return self.conduction_edge_ev - self.valence_edge_ev


@dataclass(frozen=True, slots=True)
class DopantLevel:
    """One donor or acceptor level using a declared degeneracy convention."""

    label: str
    kind: str
    level_energy_ev: float
    concentration_cm3: float
    degeneracy_factor: float
    convention: str
    source: str

    def __post_init__(self) -> None:
        if self.kind not in {"donor", "acceptor"}:
            raise ValueError("kind must be 'donor' or 'acceptor'")
        if self.concentration_cm3 < 0.0:
            raise ValueError("dopant concentration cannot be negative")
        if self.degeneracy_factor <= 0.0:
            raise ValueError("degeneracy_factor must be positive")


@dataclass(frozen=True, slots=True)
class SemiconductorResult:
    """One equilibrium charge-neutrality result with validation evidence."""

    temperature_k: float
    chemical_potential_ev: float
    electron_concentration_cm3: float
    hole_concentration_cm3: float
    ionized_donors_cm3: float = 0.0
    ionized_acceptors_cm3: float = 0.0
    neutrality_residual_cm3: float = 0.0
    metadata: dict[str, object] = field(default_factory=dict)


def fermi_dirac(
    energy_ev: ArrayLike,
    chemical_potential_ev: float,
    thermal_energy_ev: float,
) -> NDArray[np.float64]:
    """Return a stable Fermi–Dirac occupation for positive ``k_B T``."""
    if thermal_energy_ev <= 0.0:
        raise ValueError("thermal_energy_ev must be positive")
    argument = (
        np.asarray(energy_ev, dtype=np.float64) - chemical_potential_ev
    ) / thermal_energy_ev
    result = np.empty_like(argument)
    positive = argument >= 0.0
    exp_negative = np.exp(-argument[positive])
    result[positive] = exp_negative / (1.0 + exp_negative)
    exp_positive = np.exp(argument[~positive])
    result[~positive] = 1.0 / (1.0 + exp_positive)
    return result


def dos_conduction(energy_ev: ArrayLike, model: BandEdgeModel) -> NDArray[np.float64]:
    """Return the 3D conduction DOS in the documented output units."""
    raise NotImplementedError("Part II: implement unit-checked multi-valley DOS")


def dos_valence(energy_ev: ArrayLike, model: BandEdgeModel) -> NDArray[np.float64]:
    """Return the 3D valence DOS in the documented output units."""
    raise NotImplementedError("Part II: implement unit-checked multi-band DOS")


def electron_concentration(
    chemical_potential_ev: float,
    temperature_k: float,
    model: BandEdgeModel,
) -> float:
    """Return equilibrium electron concentration in cm^-3."""
    raise NotImplementedError("Part II: implement converged Fermi-Dirac integration")


def hole_concentration(
    chemical_potential_ev: float,
    temperature_k: float,
    model: BandEdgeModel,
) -> float:
    """Return equilibrium hole concentration in cm^-3."""
    raise NotImplementedError("Part II: implement converged Fermi-Dirac integration")


def solve_chemical_potential(
    temperature_k: float,
    model: BandEdgeModel,
    dopants: tuple[DopantLevel, ...] = (),
) -> SemiconductorResult:
    """Solve bracketed equilibrium charge neutrality."""
    raise NotImplementedError("Parts II/IV: implement and validate neutrality solve")
