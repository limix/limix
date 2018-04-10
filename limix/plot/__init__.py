r"""
******************
Plotting utilities
******************

This module is implemented by limix-plot_ package, from which we export
the following functions and classes:



.. _limix-plot: https://limix-plot.readthedocs.io/
"""
from __future__ import absolute_import as _
from limix_plot import (kinship, qqplot, box_aspect, image, manhattan, pca,
                        normal, power, ConsensusCurve)

__all__ = ['kinship', 'qqplot', 'box_aspect', 'image', 'manhattan', 'pca',
           'normal', 'power', 'ConsensusCurve']
