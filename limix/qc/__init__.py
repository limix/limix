r"""
***************
Quality control
***************

This module provides tools to help you normalise your data according to your
needs, such as `Box-Cox transformation`_ (:func:`limix.qc.boxcox`),
missing value imputation (:func:`limix.qc.mean_impute`), zero-mean
standardisation (:func:`limix.qc.mean_standardize`), and Normal quantile
normalisation (:func:`limix.qc.quantile_gaussianize`).


.. _Box-Cox transformation: https://en.wikipedia.org/wiki/Power_transform#Box.E2.80.93Cox_transformation

Interface
^^^^^^^^^^^^^^

.. autofunction:: limix.qc.boxcox
.. autofunction:: limix.qc.mean_impute
.. autofunction:: limix.qc.mean_standardize
.. autofunction:: limix.qc.quantile_gaussianize
.. autofunction:: limix.qc.indep_pairwise
.. autofunction:: limix.qc.compute_maf
.. autofunction:: limix.qc.count_missingness
.. autofunction:: limix.qc.regress_out
.. autofunction:: limix.qc.remove_dependent_cols
"""

from .allele import compute_maf
from .impute import mean_impute
from .ld import indep_pairwise
from .linalg import remove_dependent_cols
from .missing import count_missingness
from .regress import regress_out
from .trans import boxcox, mean_standardize, quantile_gaussianize

__all__ = [
    'boxcox', 'mean_standardize', 'quantile_gaussianize', 'regress_out',
    'remove_dependent_cols', 'mean_impute', 'indep_pairwise',
    'count_missingness', 'compute_maf'
]
