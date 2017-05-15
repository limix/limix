r"""
******************
Plotting utilities
******************

Manhattan
^^^^^^^^^

.. autofunction:: limix.plot.plot_manhattan

QQ plot
^^^^^^^

.. autofunction:: limix.plot.qqplot

Kinship
^^^^^^^

.. autofunction:: limix.plot.plot_kinship

Normal distribution
^^^^^^^^^^^^^^^^^^^

.. autofunction:: limix.plot.plot_normal

"""

from .plot_normal import plot_normal
from .qqplot import qqplot
from .manhattan import plot_manhattan
from .kinship import plot_kinship

__all__ = ['plot_manhattan', 'qqplot', 'plot_normal', 'plot_kinship']
