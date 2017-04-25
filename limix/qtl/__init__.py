r"""
**********************
Single-variant testing
**********************

- :func:`.qtl_test_lm`
- :func:`.qtl_test_lmm`
- :func:`.qtl_test_glmm`
- :func:`.qtl_test_lmm_kronecker`
- :func:`.qtl_test_interaction_lmm_kronecker`
- :func:`.qtl_test_interaction_lmm`
- :func:`.forward_lmm`
- :class:`.LMM`

Public interface
^^^^^^^^^^^^^^^^
"""

from .qtl import qtl_test_lm
from .qtl import qtl_test_lmm
from .qtl import qtl_test_lmm_kronecker
from .qtl import qtl_test_interaction_lmm_kronecker
from .qtl import qtl_test_interaction_lmm
from .qtl import forward_lmm
from .glmm import qtl_test_glmm
from .lmm import LMM

__all__ = ['qtl_test_lm', 'qtl_test_lmm', 'qtl_test_lmm_kronecker',
           'qtl_test_interaction_lmm_kronecker', 'qtl_test_interaction_lmm',
           'forward_lmm', 'LMM', 'qtl_test_glmm']
