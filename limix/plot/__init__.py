r"""
******************
Plotting utilities
******************


Box plot
^^^^^^^^

.. automethod:: limix.plot.LimixPlot.boxplot

Curve
^^^^^

.. automethod:: limix.plot.LimixPlot.curve
.. autoclass:: limix.plot.ConsensusCurve

Manhattan plot
^^^^^^^^^^^^^^

.. automethod:: limix.plot.LimixPlot.manhattan

Quantile-Quantile plot
^^^^^^^^^^^^^^^^^^^^^^

.. automethod:: limix.plot.LimixPlot.qqplot

Power plot
^^^^^^^^^^

.. automethod:: limix.plot.LimixPlot.power

Principal components plot
^^^^^^^^^^^^^^^^^^^^^^^^^

.. automethod:: limix.plot.LimixPlot.pca

Kinship plot
^^^^^^^^^^^^

.. automethod:: limix.plot.LimixPlot.kinship

Normal distribution
^^^^^^^^^^^^^^^^^^^

.. automethod:: limix.plot.LimixPlot.normal

File plot
^^^^^^^^^

.. automethod:: limix.plot.see_image

"""

from .consensus import ConsensusCurve
from .image import see_image
from .plotter import get_plot as get
from .plotter import LimixPlot

__all__ = ['ConsensusCurve', 'see_image', 'get', 'LimixPlot']
