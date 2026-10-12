# DLSU SEMPHY01

Course repository for DLSU course `SEMPHY01`.

Historical source material is preserved without normalization under `legacy/`.

## Current assessment

- [`assessments/long-assessment-1.md`](assessments/long-assessment-1.md) — broad, individually scoped Long Assessment 1 brief for the four-student course.
- [`templates/long-assessment-1-report-template.tex`](templates/long-assessment-1-report-template.tex) — editable LaTeX report skeleton.
- [`templates/long-assessment-1-report-template.pdf`](templates/long-assessment-1-report-template.pdf) — compiled report skeleton for preview.

## Assessment notebooks

- [`notebooks/01_bloch2d_validation.ipynb`](notebooks/01_bloch2d_validation.ipynb) — executed PhysKit checkerboard-band, effective-mass, 2D DOS, and intrinsic sheet-statistics validation.
- [`notebooks/02_intrinsic_effmass.ipynb`](notebooks/02_intrinsic_effmass.ipynb) — executable Part II skeleton for the 3D intrinsic model and silicon control.
- [`notebooks/03_material_analysis.ipynb`](notebooks/03_material_analysis.ipynb) — executable Part III declaration, provenance, comparison, and sensitivity skeleton.
- [`notebooks/04_doping_analysis.ipynb`](notebooks/04_doping_analysis.ipynb) — executable Part IV incomplete-ionization and doping-sweep skeleton.

Supporting scaffolds are in `src/`, `tests/`, `data/`, and `figures/`. Skeleton
notebooks and skipped test cases are explicitly marked TODO and are not yet
submission-ready.

The notebooks use Python 3.14+ and the PhysKit commit pinned in
[`requirements.txt`](requirements.txt). Create an environment with
`python3.14 -m pip install -r requirements.txt`, then run the notebook from a
clean kernel. With `uv`, launch one directly with:

```bash
uv run --python 3.14 --with-requirements requirements.txt jupyter lab
```
