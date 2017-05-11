r"""
*************
Dataset utils
*************

- :class:`.BedReader`
- :class:`.query_and`
- :class:`.build_geno_query`

Public interface
^^^^^^^^^^^^^^^^
"""

from .bed_reader import BedReader
from .query_utils import query_and
from .query_utils import build_geno_query
from .geno_iterator import GIter

__all__ = ['BedReader', 'query_and', 'build_geno_query']
