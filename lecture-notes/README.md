# SEMPHY01 Lecture Notes

Development manuscript for `dlsu_AY20262026T1`.

The manuscript develops equilibrium semiconductor carrier physics from a
solid-state and statistical-mechanical foundation. The sequence begins with the
particle in a box, Bloch states, molecular orbitals, tight-binding bands,
spectral multiplicity, density of states, partition functions, and the
chemical-potential constraint before treating intrinsic and extrinsic material.
It follows Sze and Lee's equilibrium-carrier sequence, enriched by Sze and Ng's
Fermi-integral and degeneracy discussion and Rudan's derivations from
periodic-band curvature, reciprocal-space state counting, and effective-mass
tensors.

## Student figure workflow

Each published figure has one student-facing notebook.  Notebook paths follow
this convention:

```text
figures/0X_chapter-name/fig_0X_0X.ipynb
```

The notebook named in a figure-caption footnote is the authoritative
reproduction record for that figure.  Each notebook begins with a transparent,
direct calculation from the governing equations.  This first calculation uses
NumPy explicitly so students can inspect the physical construction before any
course helper is introduced.  Its final section adds `figures/src/` to
`sys.path` and demonstrates the reusable implementation that writes the
publication PDF.

Reusable scientific and plotting functions live in:

```text
figures/src/semphy/
```

Their tests mirror the source tree module by module:

```text
figures/src/semphy/band_foundations.py
figures/tests/semphy/test_band_foundations.py
```

The same source-to-test correspondence applies to every module.  Source modules
use deliberately extensive NumPy-style docstrings.  These docstrings state the
physical model, quantities and units, normalization convention, assumptions,
parameters, returned values, exceptions, generated filename, and numerical
limitations wherever they apply.  They are written so a scientist can inspect
the physical meaning of a reusable function without first reverse-engineering
its implementation.

This separation keeps the notebooks explanatory while giving repeated
calculations a documented, tested implementation.  Students should read and
execute a figure notebook from top to bottom before adapting the reusable
function.

## Build and verify

From the repository-root `lecture-notes/` directory, execute the figure
notebooks and then compile the manuscript:

```bash
for notebook in $(find figures -mindepth 2 -name 'fig_*.ipynb' | sort); do
  uv run --python 3.14 --with jupyter --with numpy --with matplotlib \
    jupyter nbconvert --execute --to notebook --inplace "$notebook"
done

PYTHONPATH=figures/src \
  uv run --python 3.14 --with 'pytest>=8,<10' --with numpy \
  pytest -q figures/tests

latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The manuscript output is `main.pdf`; generated figure PDFs remain directly
under `figures/`.

## Current chapters

1. **From a Particle in a Box to Bloch States** — hard-wall confinement,
   spectral representation, periodic boundary conditions, translation
   symmetry, Bloch's theorem, reciprocal vectors, and the Brillouin zone.
2. **From HOMO--LUMO Levels to Tight-Binding Bands** — localized orbitals,
   bonding and antibonding states, finite-chain molecular levels, the Bloch
   basis, tight-binding dispersion, multiorbital bands, and gap formation.
3. **Spectral Representation and the Density of States** — spectral
   projectors, discrete multiplicity, Bloch bands, reciprocal-space state
   counting, the general band DOS, constant-energy surfaces, van Hove
   singularities, dimensionality, anisotropic DOS masses, numerical quadrature,
   and normalization.
4. **Statistical Mechanics of Carrier Populations** — equipartition and its
   limits, reservoir multiplicity, canonical and grand partition functions,
   one-electron band probabilities, many-electron mean occupations, the
   Boltzmann limit, and carrier integrals.
5. **Determining the Chemical Potential from the Available States** — fixed
   particle number, charge neutrality, energy alignment, model selection,
   monotonic closure equations, spatially varying bands, and the required
   physical inputs.
6. **The Intrinsic Semiconductor** — band formation, full-band and hole
   interpretation, band-edge expansions, tensor effective masses,
   Fermi-integral populations, the nondegenerate limit, intrinsic chemical
   potential, intrinsic concentration, mass action, temperature dependence,
   and a silicon control estimate.
7. **The Extrinsic Semiconductor** — substitutional donors and acceptors,
   hydrogenic binding, localized-state degeneracy, incomplete ionization,
   charge neutrality, compensation, majority and minority carriers, freeze-out,
   temperature regimes, degeneracy, band-gap narrowing, nonuniform doping, and
   a verified computational neutrality solve.

Bibliography metadata is in `bibliography/references.bib`.  Every figure
caption carries a footnote to its corresponding notebook, allowing a student
to move directly from the rendered physical claim to the calculation that
produced it.
