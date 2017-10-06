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
from .extract import extract
from .set_utils import annotate_sets, sets_from_bim
from .temp import TemporaryDirectory
from .timer import Timer
from .url import download, urlretrieve

log2pi = 1.837877066409345339081937709124758839607238769531250

__all__ = [
    'sets_from_bim', 'annotate_sets', 'TemporaryDirectory', 'urlretrieve',
    'Timer', 'npy_dask', 'array_hash', 'log2pi', 'download', 'extract',
    'filehash'
]
