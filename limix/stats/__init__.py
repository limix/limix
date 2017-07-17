r"""
**********
Statistics
**********

.. autofunction:: limix.stats.effsizes_se
.. autofunction:: limix.stats.lrt_pvalues

PCA
^^^

.. autofunction:: limix.stats.pca


P-values
^^^^^^^^

.. autofunction:: limix.stats.qvalues
.. autofunction:: limix.stats.empirical_pvalues

Kinship processing
^^^^^^^^^^^^^^^^^^

.. autofunction:: limix.stats.linear_kinship
.. autofunction:: limix.stats.gower_norm
.. autofunction:: limix.stats.indep_pairwise
.. autofunction:: limix.stats.maf

Chi2
^^^^

.. autofunction:: limix.stats.Chi2mixture

"""

from .chi2mixture import Chi2mixture
from .confusion import confusion_matrix
from .fdr import qvalues
from .kinship import gower_norm, linear_kinship
from .lrt import effsizes_se, lrt_pvalues
from .pca import pca
from .preprocess import compute_maf, count_missingness, indep_pairwise
from .teststats import empirical_pvalues

__all__ = [
    'pca', 'boxcox', 'gower_norm', 'qvalues', 'empirical_pvalues',
    'Chi2mixture', 'indep_pairwise', 'maf', 'linear_kinship', 'lrt_pvalues',
    'effsizes_se', 'count_missingness', 'confusion_matrix'
]
