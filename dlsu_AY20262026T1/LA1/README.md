# SEMPHY01 Long Assessment 1

This directory is the self-contained Long Assessment 1 workspace for
`dlsu_AY20262026T1`.

## Assessment brief

- [`long-assessment-1.md`](long-assessment-1.md) — assessment scope,
  requirements, evidence expectations, and submission guidance.

## Computational notebooks

- [`notebooks/01_bloch2d_validation.ipynb`](notebooks/01_bloch2d_validation.ipynb)
  — executed checkerboard-band, effective-mass, two-dimensional DOS, and
  intrinsic sheet-statistics validation.
- [`notebooks/02_intrinsic_effmass.ipynb`](notebooks/02_intrinsic_effmass.ipynb)
  — executable Part II scaffold for the three-dimensional intrinsic model and
  silicon control.
- [`notebooks/03_material_analysis.ipynb`](notebooks/03_material_analysis.ipynb)
  — Part III declaration, provenance, comparison, and sensitivity scaffold.
- [`notebooks/04_doping_analysis.ipynb`](notebooks/04_doping_analysis.ipynb)
  — Part IV incomplete-ionization and doping-sweep scaffold.

Skeleton notebooks and skipped test cases are explicitly marked and are not yet
submission-ready.

## Supporting files

- `src/` — reusable equilibrium-semiconductor model scaffold.
- `tests/` — validation tests and explicitly skipped implementation targets.
- `data/` — material parameters, source notes, and digitized band-data scaffold.
- `figures/` — generated assessment figures.
- `templates/` — editable and compiled LaTeX report skeleton.
- `requirements.txt` — Python dependencies and pinned PhysKit commit.
- `pytest.ini` — test discovery and local import configuration.

## Python environment and tests

Run commands from this directory:

```bash
python3.14 -m pip install -r requirements.txt
pytest
```

With `uv`:

```bash
uv run --python 3.14 --with-requirements requirements.txt pytest
uv run --python 3.14 --with-requirements requirements.txt jupyter lab
```

## Report template

Build the report from its own directory so bibliography and figure paths resolve
consistently:

```bash
cd templates
latexmk -pdf long-assessment-1-report-template.tex
```

The compiled PDF is tracked so students can inspect the skeleton without a
local TeX installation.
