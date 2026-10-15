# DLSU SEMPHY01

Course repository for DLSU course `SEMPHY01`.

## `dlsu_AY20262026T1`

Long Assessment 1 is consolidated in
[`dlsu_AY20262026T1/LA1/`](dlsu_AY20262026T1/LA1/). Its assessment brief,
notebooks, source scaffold, tests, data, figures, dependency pin, and report
template are all contained under that directory.

Start from the assessment root so notebook imports and relative paths remain
consistent:

```bash
cd dlsu_AY20262026T1/LA1
python3.14 -m pip install -r requirements.txt
pytest
```

With `uv`, launch JupyterLab from the same directory:

```bash
cd dlsu_AY20262026T1/LA1
uv run --python 3.14 --with-requirements requirements.txt jupyter lab
```

See the [LA1 workspace guide](dlsu_AY20262026T1/LA1/README.md) for the complete
layout and build commands.

The developing course manuscript is at the repository root in
[`lecture-notes/`](lecture-notes/). Its opening sequence connects the particle
in a box to Bloch states and HOMO--LUMO levels to tight-binding bands before
introducing density of states and equilibrium carrier statistics.  Each
published figure has a student-facing reproduction notebook identified by a
footnote in its caption; see the
[lecture-note workflow](lecture-notes/README.md).

## Historical material

Historical source material is preserved without normalization under
[`legacy/`](legacy/).
