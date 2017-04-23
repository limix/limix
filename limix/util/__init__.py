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

from .set_utils import sets_from_bim
from .set_utils import annotate_sets
from .geno_utils import estCumPos
from .geno_utils import unique_variants
from .temp import TemporaryDirectory
from .url import urlretrieve
from . import preprocess

__all__ = ['sets_from_bim',
           'annotate_sets',
           'estCumPos',
           'unique_variants',
           'TemporaryDirectory',
           'urlretrieve']
