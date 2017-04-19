r"""
**********************
Single-variant testing
**********************

- :func:`.qtl_test_lm`
- :func:`.qtl_test_lmm`
- :class:`.LMM`

Public interface
^^^^^^^^^^^^^^^^
"""

from .qtl import qtl_test_lm
from .qtl import qtl_test_lmm
from .lmm import LMM

__all__ = ['qtl_test_lm', 'qtl_test_lmm', 'LMM']

