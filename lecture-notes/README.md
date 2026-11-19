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
8. **Brownian Motion, Diffusion, and the Heat Equation** — random walks,
   Wiener processes, Itô calculus, Fokker--Planck evolution, Fick's law,
   finite-difference Laplacians, matrix-exponential propagation, and the heat
   equation.
9. **From Langevin Dynamics to the Maxwell--Boltzmann Gas** — inertial thermal
   motion, fluctuation--dissipation balance, Maxwell velocities, velocity
   memory, and the ballistic-to-Brownian crossover.
10. **The Drude and Drude--Sommerfeld Electron Gases** — DC and AC Drude
    response, mobility, the Einstein relation, Fermi--Dirac occupation, and the
    Fermi-window conductivity tensor.
11. **Effective Mass and the Semiconductor Carrier Model** — semiclassical
    band dynamics, effective-mass tensors, multivalley transport, electron and
    hole drift--diffusion currents, and carrier continuity.
12. **Mobility, Resistivity, and the Hall Effect** — lattice and impurity
    scattering, Matthiessen's rule, conductivity, resistivity, Hall carrier
    sign, Hall density, and Hall mobility.
13. **Generation, Recombination, and Carrier Continuity** — excess-carrier
    lifetime, direct and trap-assisted recombination, quasi-Fermi levels,
    diffusion length, surface recombination, and Haynes--Shockley transport.
14. **Barrier Transport and Space-Charge-Limited Current** — thermionic
    emission, rectangular-barrier tunneling, injected space charge, and the
    constant-mobility and velocity-saturated current laws.
15. **High-Field Carrier Transport** — velocity saturation, transferred
    electrons, negative differential mobility, and impact ionization.
16. **Equilibrium p--n Junction and Depletion Electrostatics** — contact
    potential, charge neutrality, abrupt and graded depletion layers, and a
    nonlinear Poisson--Boltzmann validation.
17. **Depletion Capacitance and Doping Profiling** — incremental charge,
    abrupt and graded C--V laws, varactors, and numerical profile extraction.
18. **Junction Current--Voltage Characteristics** — the law of the junction,
    minority diffusion, the ideal diode equation, depletion-region
    generation--recombination, high injection, and series resistance.
19. **Charge Storage, Diffusion Capacitance, and Switching** — charge control,
    reverse recovery, and implicit transient minority diffusion.
20. **Junction Breakdown and Heterojunctions** — tunneling, avalanche
    integration, field crowding, band offsets, and heterojunction voltage
    partition.

## Sze--Lee Chapter 2 computational sequence

Four executed notebooks complete the semiconductor transport sequence:

- [`fig_12_01.ipynb`](figures/12_mobility-resistivity-hall/fig_12_01.ipynb)
  connects scattering-limited mobility to resistivity and Hall response.
- [`fig_13_01.ipynb`](figures/13_generation-recombination-continuity/fig_13_01.ipynb)
  verifies lifetime decay, quasi-Fermi splitting, diffusion length, and moving
  minority-carrier packets.
- [`fig_14_01.ipynb`](figures/14_barrier-space-charge-transport/fig_14_01.ipynb)
  compares thermionic activation, tunneling, and space-charge current laws.
- [`fig_15_01.ipynb`](figures/15_high-field-transport/fig_15_01.ipynb)
  evaluates velocity saturation, intervalley transfer, and impact ionization.

## Sze--Lee Chapter 3 computational sequence

Five executed notebooks develop the p--n junction:

- [`fig_16_01.ipynb`](figures/16_equilibrium-pn-junction/fig_16_01.ipynb)
  compares depletion electrostatics with nonlinear Poisson--Boltzmann finite
  differences.
- [`fig_17_01.ipynb`](figures/17_depletion-capacitance-profiling/fig_17_01.ipynb)
  validates C--V linearization and one-sided doping extraction.
- [`fig_18_01.ipynb`](figures/18_junction-current-voltage/fig_18_01.ipynb)
  checks finite-base minority diffusion and its junction-edge current.
- [`fig_19_01.ipynb`](figures/19_charge-storage-switching/fig_19_01.ipynb)
  resolves transient minority storage and extraction with implicit finite
  differences.
- [`fig_20_01.ipynb`](figures/20_breakdown-heterojunctions/fig_20_01.ipynb)
  computes avalanche breakdown and ideal heterojunction band offsets.

The notebooks also generate thirteen original schematic plates covering every
physical figure subject in Sze--Lee Chapter 3, Figs. 1--30: junction formation,
band diagrams, depletion profiles, bias response, capacitance construction,
injection, diode characteristics, temperature response, the small-signal
circuit, switching, breakdown, curvature, and heterojunction lineup. The
schematics are generated from the chapter equations and do not reproduce
scanned textbook artwork.

## Brownian-motion transport labs

Five executed, self-contained computational labs connect microscopic random
motion to macroscopic transport:

- [`lab-01-brownian-to-ficks-law.ipynb`](figures/08_brownian-motion-diffusion/lab-01-brownian-to-ficks-law.ipynb)
  develops mean-square displacement, the ensemble density, and Fick's laws.
- [`lab-02-particle-current.ipynb`](figures/08_brownian-motion-diffusion/lab-02-particle-current.ipynb)
  measures signed particle crossings and verifies Fick current and local
  continuity.
- [`lab-03-matrix-diffusion-exponential-evolution.ipynb`](figures/08_brownian-motion-diffusion/lab-03-matrix-diffusion-exponential-evolution.ipynb)
  constructs the heat/diffusion matrix from discrete fluxes and evolves the
  system with a matrix exponential.
- [`lab-04-langevin-maxwell-boltzmann-drude.ipynb`](figures/08_brownian-motion-diffusion/lab-04-langevin-maxwell-boltzmann-drude.ipynb)
  connects inertial Brownian motion to Maxwell–Boltzmann equilibrium, the
  Einstein relation, and DC and AC Drude response.
- [`lab-05-effective-mass-from-bands-to-transport.ipynb`](figures/08_brownian-motion-diffusion/lab-05-effective-mass-from-bands-to-transport.ipynb)
  extracts the effective-mass tensor from band curvature and distinguishes its
  density-of-states and multivalley transport roles.

Bibliography metadata is in `bibliography/references.bib`.  Every figure
caption carries a footnote to its corresponding notebook, allowing a student
to move directly from the rendered physical claim to the calculation that
produced it.
