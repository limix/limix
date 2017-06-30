r"""
**********************
Single-variant testing
**********************


Linear mixed models
^^^^^^^^^^^^^^^^^^^

.. autofunction:: limix.lmm.qtl_test_lmm

Generalised linear mixed models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: limix.glmm.qtl_test_glmm

"""

from .glmm import qtl_test_glmm
from .lmm import qtl_test_lmm

__all__ = ['qtl_test_lmm', 'qtl_test_glmm']
