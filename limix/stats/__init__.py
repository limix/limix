r"""
**********
Statistics
**********

- :func:`.pca`
- :func:`.boxcox`

Public interface
^^^^^^^^^^^^^^^^
"""

from .pca import pca
from .trans import boxcox

__all__ = ['pca', 'boxcox']
