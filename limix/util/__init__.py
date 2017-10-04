r"""
*****************
Utility functions
*****************

- :func:`.sets_from_bim`
- :func:`.annotate_sets`
- :func:`.estCumPos`
- :func:`.unique_variants`
- :class:`.TemporaryDirectory`
- :func:`.urlretrieve`

Public interface
^^^^^^^^^^^^^^^^
"""

from . import npy_dask
from ._hash import array_hash
from .geno_utils import estCumPos, unique_variants
from .set_utils import annotate_sets, sets_from_bim
from .temp import TemporaryDirectory
from .timer import Timer
from .url import urlretrieve

log2pi = 1.837877066409345339081937709124758839607238769531250

__all__ = [
    'sets_from_bim', 'annotate_sets', 'estCumPos', 'unique_variants',
    'TemporaryDirectory', 'urlretrieve', 'Timer', 'npy_dask', 'array_hash',
    'log2pi'
]
