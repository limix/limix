r"""
*************
Dataset utils
*************

- :class:`.BedReader`
- :class:`.QueryList`

Public interface
^^^^^^^^^^^^^^^^
"""

from .bed_reader import BedReader
from .query import QueryList
from .geno_iterator import GIter

__all__ = ['BedReader', 'QueryList']
