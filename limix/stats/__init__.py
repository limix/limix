r"""
**********
Statistics
**********

- :func:`.pca`
- :func:`.boxcox`
- :func:`.gower_norm`

Public interface
^^^^^^^^^^^^^^^^
"""

from .pca import pca
from .trans import boxcox
from .kinship import gower_norm

__all__ = ['pca', 'boxcox', 'gower_norm']
