"""
Statistics.
"""

from ._allele import allele_expectation, allele_frequency, compute_dosage
from ._chi2 import Chi2Mixture
from ._confusion import confusion_matrix
from ._kinship import linear_kinship
from ._lrt import lrt_pvalues
from ._pca import pca
from ._pvalue import empirical_pvalues, multipletests

__all__ = [
    "Chi2Mixture",
    "allele_expectation",
    "allele_frequency",
    "compute_dosage",
    "confusion_matrix",
    "empirical_pvalues",
    "linear_kinship",
    "lrt_pvalues",
    "multipletests",
    "pca",
]
