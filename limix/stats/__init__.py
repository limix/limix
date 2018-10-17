r"""
**********
Statistics
**********

Principal component analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: limix.stats.pca
    :noindex:

P-value correction
^^^^^^^^^^^^^^^^^^

.. autofunction:: limix.stats.multipletests
    :noindex:
.. autofunction:: limix.stats.empirical_pvalues
    :noindex:
.. autofunction:: limix.stats.Chi2Mixture
    :noindex:

Ground truth evaluation
^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: limix.stats.confusion_matrix
    :noindex:

Kinship
^^^^^^^

.. autofunction:: limix.stats.linear_kinship
    :noindex:

Likelihood-ratio test
^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: limix.stats.effsizes_se
    :noindex:
.. autofunction:: limix.stats.lrt_pvalues
    :noindex:
"""

from ._chi2mixture import Chi2Mixture
from ._confusion import confusion_matrix
from ._kinship import linear_kinship
from ._lrt import effsizes_se, lrt_pvalues
from ._allele import allele_frequency, compute_dosage, allele_expectation
from ._pvalue import multipletests, empirical_pvalues
from ._pca import pca

__all__ = [
    "pca",
    "multipletests",
    "empirical_pvalues",
    "Chi2Mixture",
    "linear_kinship",
    "lrt_pvalues",
    "effsizes_se",
    "confusion_matrix",
    "allele_frequency",
    "compute_dosage",
    "allele_expectation",
]
