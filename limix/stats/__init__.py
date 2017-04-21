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

Public interface
^^^^^^^^^^^^^^^^
"""

from .pca import pca
from .trans import boxcox
from .kinship import gower_norm
from .fdr import qvalues
from .teststats import empirical_pvalues
from .chi2mixture import Chi2mixture

__all__ = ['pca', 'boxcox', 'gower_norm', 'qvalues',
           'empirical_pvalues', 'Chi2mixture']
