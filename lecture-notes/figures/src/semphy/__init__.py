"""Scientific calculations and reproducible figures for SEMPHY01.

The package is intentionally small and readable.  Student notebooks first
construct each result directly from the governing equations, then add
``figures/src`` to ``sys.path`` and import a reusable implementation from this
package.  This sequence keeps the physics visible before introducing software
abstraction.

Modules
-------
band_foundations
    Fixed figure recipes for confinement, Bloch states, molecular levels, and
    tight-binding bands.
density_of_states
    Unit-aware numerical relations for sampled bands, analytic densities of
    states, effective masses, broadening, and normalization.
density_of_states_figures
    Fixed figure recipes that turn the DOS relations into Chapter 3 graphics.
equilibrium_semiconductors
    Fixed figure recipes for intrinsic populations, dopants, compensation, and
    temperature-dependent charge neutrality.

The mirrored test tree under ``figures/tests/semphy`` gives every source module
an immediately discoverable test file.  These implementations are educational
and inspectable; they do not claim the convergence controls or material models
required of production electronic-structure or device-simulation software.
"""

from __future__ import annotations
