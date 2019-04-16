"""
Statistics.
"""

from ._chi2 import Chi2Mixture
from ._confusion import confusion_matrix
from ._kinship import linear_kinship
from ._lrt import effsizes_se, lrt_pvalues
from ._allele import allele_frequency, compute_dosage, allele_expectation
from ._pvalue import multipletests, empirical_pvalues
from ._pca import pca

__all__ = [
    "Chi2Mixture",
    "allele_expectation",
    "allele_frequency",
    "compute_dosage",
    "confusion_matrix",
    "effsizes_se",
    "empirical_pvalues",
    "linear_kinship",
    "lrt_pvalues",
    "multipletests",
    "pca",
]
