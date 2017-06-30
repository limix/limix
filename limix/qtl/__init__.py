r"""
**********************
Single-variant testing
**********************

Linear models
^^^^^^^^^^^^^

.. autofunction:: limix.qtl.qtl_test_lm

Linear mixed models
^^^^^^^^^^^^^^^^^^^

.. autofunction:: limix.lmm.qtl_test_lmm
.. autofunction:: limix.qtl.qtl_test_interaction_lmm
.. autofunction:: limix.qtl.qtl_test_lmm_kronecker
.. autofunction:: limix.qtl.qtl_test_interaction_lmm_kronecker
.. autofunction:: limix.qtl.forward_lmm

Generalised linear mixed models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: limix.glmm.qtl_test_glmm

"""

from .glmm import qtl_test_glmm
from .lmm import qtl_test_lmm
from .qtl import (forward_lmm, qtl_test_interaction_lmm,
                  qtl_test_interaction_lmm_kronecker, qtl_test_lm,
                  qtl_test_lmm_kronecker)

__all__ = [
    'qtl_test_lm', 'qtl_test_lmm', 'qtl_test_lmm_kronecker',
    'qtl_test_interaction_lmm_kronecker', 'qtl_test_interaction_lmm',
    'forward_lmm', 'qtl_test_glmm'
]
