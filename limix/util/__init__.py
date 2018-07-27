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
from ._hash import array_hash, filehash
from ._session_text import session_text
from .extract import extract
from .set_utils import annotate_sets, sets_from_bim
from .timer import Timer
from .tmpdir import TmpDir
from .url import download, urlretrieve

log2pi = 1.837877066409345339081937709124758839607238769531250

__all__ = [
    "annotate_sets",
    "array_hash",
    "download",
    "extract",
    "filehash",
    "log2pi",
    "npy_dask",
    "session_text",
    "sets_from_bim",
    "Timer",
    "TmpDir",
    "urlretrieve",
]
