r"""Plotting utilities"""

from .consensus import ConsensusCurve
from .limix_plot import LimixPlot

plot = LimixPlot()

__all__ = ['LimixPlot', 'plot', 'ConsensusCurve']
