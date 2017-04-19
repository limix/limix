r"""
**********
Statistics
**********

- :func:`.pca`
- :func:`.boxcox`
- :func:`.gower_norm`
- :class:`.Chi2mixture`

Public interface
^^^^^^^^^^^^^^^^
"""

from .pca import pca
from .trans import boxcox
from .kinship import gower_norm
from .chi2mixture import Chi2mixture

__all__ = ['pca', 'boxcox', 'gower_norm', 'Chi2mixture']
