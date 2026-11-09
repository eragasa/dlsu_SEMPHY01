"""Compatibility import for the notebook-colocated semiconductor model.

The canonical implementation is ``notebooks/semiconductor_model.py`` so every
notebook and its required local library travel together.  Existing tests retain
the historical ``src.semiconductor_model`` import through this re-export.
"""

from notebooks.semiconductor_model import *  # noqa: F403
