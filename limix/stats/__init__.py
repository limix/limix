r"""
**********
Statistics
**********

Principal component analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: limix.stats.pca

P-value correction
^^^^^^^^^^^^^^^^^^

.. autofunction:: limix.stats.multipletests
.. autofunction:: limix.stats.empirical_pvalues
.. autofunction:: limix.stats.Chi2Mixture

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

from .chi2mixture import Chi2Mixture
from ._confusion import confusion_matrix
from .kinship import linear_kinship
from .lrt import effsizes_se, lrt_pvalues
from .pca import pca
from .pvalue import multipletests
from .teststats import empirical_pvalues

__all__ = [
    "pca",
    "multipletests",
    "empirical_pvalues",
    "Chi2Mixture",
    "linear_kinship",
    "lrt_pvalues",
    "effsizes_se",
    "confusion_matrix",
]
