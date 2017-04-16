r"""
************************
Quantitative trait locus
************************

- :func:`.forward_lmm`
- :func:`.forward_lmm_kronecker`
- :func:`.qtl_test_interaction_lmm`
- :func:`.qtl_test_interaction_lmm_kronecker`
- :func:`.qtl_test_lm`
- :func:`.qtl_test_lmm`
- :func:`.qtl_test_lmm_kronecker`

Public interface
^^^^^^^^^^^^^^^^
"""

from .qtl import (forward_lmm, forward_lmm_kronecker, qtl_test_interaction_lmm,
                  qtl_test_interaction_lmm_kronecker, qtl_test_lm,
                  qtl_test_lmm, qtl_test_lmm_kronecker)

__all__ = [
    'forward_lmm', 'forward_lmm_kronecker', 'qtl_test_interaction_lmm',
    'qtl_test_interaction_lmm_kronecker', 'qtl_test_lm', 'qtl_test_lmm',
    'qtl_test_lmm_kronecker'
]
