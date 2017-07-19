r"""
**********
Statistics
**********

.. autofunction:: limix.stats.effsizes_se
.. autofunction:: limix.stats.lrt_pvalues
.. autofunction:: limix.stats.confusion_matrix

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
.. autofunction:: limix.stats.compute_maf

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
from .teststats import empirical_pvalues

__all__ = [
    'pca', 'boxcox', 'gower_norm', 'qvalues', 'empirical_pvalues',
    'Chi2mixture', 'linear_kinship', 'lrt_pvalues', 'effsizes_se',
    'confusion_matrix'
]
