r"""
*****************
Utility functions
*****************

- :func:`.sets_from_bim`
- :func:`.annotate_sets`
- :class:`.TemporaryDirectory`
- :func:`.urlretrieve`

Public interface
^^^^^^^^^^^^^^^^
"""

from . import npy_dask
from .set_utils import annotate_sets, sets_from_bim
from .timer import Timer
from .tmpdir import TmpDir
from ._array_hash import array_hash

log2pi = 1.837877066409345339081937709124758839607238769531250

__all__ = [
    "annotate_sets",
    "log2pi",
    "npy_dask",
    "sets_from_bim",
    "Timer",
    "TmpDir",
    "array_hash",
]
