r"""
**********
Statistics
**********

Principal component analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: limix.stats.pca
.. autofunction:: limix.stats.empirical_pvalues
.. autofunction:: limix.stats.Chi2mixture

P-value correction
^^^^^^^^^^^^^^^^^^

.. autofunction:: limix.stats.multipletests

Ground truth evaluation
^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: limix.stats.confusion_matrix

Kinship
^^^^^^^

.. autofunction:: limix.stats.linear_kinship

Likelihood-ratio test
^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: limix.stats.effsizes_se
.. autofunction:: limix.stats.lrt_pvalues

"""

from .chi2mixture import Chi2mixture
from .confusion import confusion_matrix
from .kinship import linear_kinship
from .lrt import effsizes_se, lrt_pvalues
from .pca import pca
from .pvalue import multipletests
from .teststats import empirical_pvalues

__all__ = [
    'pca', 'multipletests', 'empirical_pvalues', 'Chi2mixture',
    'linear_kinship', 'lrt_pvalues', 'effsizes_se', 'confusion_matrix'
]
