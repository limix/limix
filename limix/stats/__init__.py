r"""
**********
Statistics
**********

PCA
^^^

.. autofunction:: limix.stats.pca

Normalization
^^^^^^^^^^^^^

.. autofunction:: limix.stats.boxcox

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
from .fdr import qvalues
from .kinship import gower_norm, linear_kinship
from .pca import pca
from .preprocess import indep_pairwise, maf
from .teststats import empirical_pvalues
from .trans import boxcox

__all__ = [
    'pca', 'boxcox', 'gower_norm', 'qvalues', 'empirical_pvalues',
    'Chi2mixture', 'indep_pairwise', 'maf', 'linear_kinship'
]
