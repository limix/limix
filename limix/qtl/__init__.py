r"""
**********************
Single-variant testing
**********************

- :func:`.qtl_test_lm`
- :class:`.LMM`

Public interface
^^^^^^^^^^^^^^^^
"""

from .qtl import qtl_test_lm
from .lmm import LMM

__all__ = ['qtl_test_lm', 'LMM']

