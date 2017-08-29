r"""
***************
Quality control
***************

Introduction
^^^^^^^^^^^^

This module provides the following functions to help normalise your data
according to your needs:

- `Box-Cox transformation`_ (:func:`limix.qc.boxcox`)
- Missing value imputation (:func:`limix.qc.mean_impute`)
- Gower rescaling for covariance matrices (:func:`limix.qc.gower_norm`)
- zero-mean Standardisation (:func:`limix.qc.mean_standardize`)
- Normal quantile normalisation (:func:`limix.qc.quantile_gaussianize`)
- Pairwise variant independence (:func:`limix.qc.indep_pairwise`)
- `Minor allele frequency`_ (:func:`limix.qc.compute_maf`)
- Number of missing values (:func:`limix.qc.count_missingness`)
- Plain linear regression (:func:`limix.qc.regress_out`)
- Removal of dependent columns (:func:`limix.qc.remove_dependent_cols`)

.. _Box-Cox transformation: https://en.wikipedia.org/wiki/Power_transform#Box.E2.80.93Cox_transformation
.. _Minor allele frequency: https://en.wikipedia.org/wiki/Minor_allele_frequency

Interface
^^^^^^^^^

.. autofunction:: limix.qc.boxcox
.. autofunction:: limix.qc.mean_impute
.. autofunction:: limix.qc.gower_norm
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
from .kinship import gower_norm
from .ld import indep_pairwise
from .linalg import remove_dependent_cols
from .missing import count_missingness
from .regress import regress_out
from .trans import boxcox, mean_standardize, quantile_gaussianize
from .unique_variants import unique_variants

__all__ = [
    'boxcox', 'mean_standardize', 'quantile_gaussianize', 'regress_out',
    'remove_dependent_cols', 'mean_impute', 'indep_pairwise',
    'count_missingness', 'compute_maf', 'gower_norm', 'unique_variants'
]
