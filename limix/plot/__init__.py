r"""
******************
Plotting utilities
******************

- :func:`.plot_manhattan`
- :func:`.plot_normal`
- :func:`.qqplot`

Public interface
^^^^^^^^^^^^^^^^
"""

from .plot_normal import plot_normal
from .qqplot import qqplot
from .manhattan import plot_manhattan

__all__ = ['plot_normal', 'qqplot', 'plot_manhattan']
