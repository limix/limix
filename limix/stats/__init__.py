r"""
**********
Statistics
**********

- :func:`.pca`
- :func:`.boxcox`
- :func:`.gower_norm`
- :func:`.qvalues`
- :class:`.Chi2mixture`

Public interface
^^^^^^^^^^^^^^^^
"""

from .pca import pca
from .trans import boxcox
from .kinship import gower_norm
from .fdr import qvalues
from .chi2mixture import Chi2mixture

__all__ = ['pca', 'boxcox', 'gower_norm', 'qvalues', 'Chi2mixture']
