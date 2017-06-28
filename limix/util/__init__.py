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

from ._array import asarray
from .geno_utils import estCumPos, unique_variants
from .set_utils import annotate_sets, sets_from_bim
from .temp import TemporaryDirectory
from .url import urlretrieve

__all__ = [
    'sets_from_bim', 'annotate_sets', 'estCumPos', 'unique_variants',
    'TemporaryDirectory', 'urlretrieve', 'asarray'
]
