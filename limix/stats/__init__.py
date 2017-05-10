r"""
**********
Statistics
**********

- :func:`.pca`
- :func:`.boxcox`
- :func:`.gower_norm`
- :func:`.qvalues`
- :func:`.empirical_pvalues`
- :class:`.Chi2mixture`
- :func:`.indep_pairwise`
- :func:`.maf`

Public interface
^^^^^^^^^^^^^^^^
"""

from .chi2mixture import Chi2mixture
from .fdr import qvalues
from .kinship import gower_norm
from .pca import pca
from .preprocess import indep_pairwise, maf
from .teststats import empirical_pvalues
from .trans import boxcox

__all__ = [
    'pca', 'boxcox', 'gower_norm', 'qvalues', 'empirical_pvalues',
    'Chi2mixture', 'indep_pairwise', 'maf'
]
