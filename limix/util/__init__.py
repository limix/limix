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

from .geno_utils import estCumPos, unique_variants
from .set_utils import annotate_sets, sets_from_bim
from .temp import TemporaryDirectory
from .url import urlretrieve
from .timer import Timer
from ._hash import array_hash
from . import npy_dask

__all__ = [
    'sets_from_bim', 'annotate_sets', 'estCumPos', 'unique_variants',
    'TemporaryDirectory', 'urlretrieve', 'Timer', 'npy_dask', 'array_hash'
]
