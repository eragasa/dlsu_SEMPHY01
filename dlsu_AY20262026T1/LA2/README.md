# SEMPHY01 Long Assessment 2

This directory is the self-contained Long Assessment 2 workspace for
`dlsu_AY20262026T1`.

## Assessment brief

- [`long-assessment-2.md`](long-assessment-2.md) — assessment scope,
  LA1-to-LA2 transfer requirements, transport and recombination models,
  p–n-junction prediction, validation, submission, and oral defense.

## Executed example notebooks

The examples use silicon and gallium arsenide to demonstrate the computational
workflow. They are reference calculations, not parameter sources for a
student's selected LA1 material and not submission templates to copy unchanged.

- [`01_si_gaas_transport.ipynb`](notebooks/01_si_gaas_transport.ipynb) —
  equilibrium-state transfer, doping-dependent mobility, conductivity,
  resistivity, Hall response, and the Einstein relation.
- [`02_si_gaas_recombination.ipynb`](notebooks/02_si_gaas_recombination.ipynb) —
  lifetime assumptions, diffusion length, steady finite-difference minority
  profiles, and transient charge removal.
- [`03_si_gaas_pn_electrostatics.ipynb`](notebooks/03_si_gaas_pn_electrostatics.ipynb) —
  built-in potential, depletion width, charge, field, potential, capacitance,
  and C–V doping extraction.
- [`04_si_gaas_diode_prediction.ipynb`](notebooks/04_si_gaas_diode_prediction.ipynb) —
  finite-base saturation current, ideal diode response, stored charge,
  diffusion capacitance, and material comparison.

Each notebook begins with direct NumPy/SciPy calculations before importing the
reusable implementation in [`notebooks/la2_models.py`](notebooks/la2_models.py).

## Data templates

- [`data/la1-equilibrium-states-template.csv`](data/la1-equilibrium-states-template.csv)
  — required LA1 state-transfer columns.
- [`data/example-si-gaas-equilibrium-states.csv`](data/example-si-gaas-equilibrium-states.csv)
  — illustrative machine-readable transfer table consumed by the examples.
- [`data/transport-parameters.yaml`](data/transport-parameters.yaml) — mobility
  and Hall parameter template.
- [`data/recombination-parameters.yaml`](data/recombination-parameters.yaml) —
  lifetime and recombination parameter template.
- [`data/device-design.yaml`](data/device-design.yaml) — junction geometry and
  operating-condition template.
- [`data/source-notes.md`](data/source-notes.md) — provenance notes template.

The silicon and GaAs values used by the example notebooks are intentionally
stored in the readable local module rather than presented as universal values.
Every student must replace them with sources appropriate to the exact material,
phase, temperature, concentration, direction, and sample condition inherited
from LA1.

## Report template

- [`templates/long-assessment-2-report-template.tex`](templates/long-assessment-2-report-template.tex)
  provides a compiling structure organized around transfer integrity, transport,
  continuity, electrostatics, current, uncertainty, and validation.
- [`templates/long-assessment-2-report-template.pdf`](templates/long-assessment-2-report-template.pdf)
  is the compiled reference copy.

Build from the template directory with `latexmk -pdf`.

## Environment

From the LA2 directory:

```bash
uv run --python 3.14 --with-requirements requirements.txt pytest
uv run --python 3.14 --with-requirements requirements.txt jupyter lab
```

The workspace must remain independent of Project Koios and PhysKit.
