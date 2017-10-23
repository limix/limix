r"""Plotting utilities"""

from .limix_plot import LimixPlot
from .consensus import ConsensusCurve

plot = LimixPlot()

__all__ = ['LimixPlot', 'plot', 'ConsensusCurve']
