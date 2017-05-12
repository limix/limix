r"""
******************
Plotting utilities
******************

- :func:`.plot_manhattan`
- :func:`.qqplot`
- :func:`.plot_normal`

Public interface
^^^^^^^^^^^^^^^^
"""

from .plot_normal import plot_normal
from .qqplot import qqplot
from .manhattan import plot_manhattan
from .kinship import plot_kinship

__all__ = ['plot_manhattan', 'qqplot', 'plot_normal', 'plot_kinship']
