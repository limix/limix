r"""
******************
Plotting utilities
******************

This module is implemented by limix-plot_ package, from which we export
the following functions and classes:

- :py:class:`~limix_plot.ConsensusCurve`
- :py:func:`~limix_plot.box_aspect`
- :py:func:`~limix_plot.image`
- :py:func:`~limix_plot.kinship`
- :py:func:`~limix_plot.load_dataset`
- :py:func:`~limix_plot.manhattan`
- :py:func:`~limix_plot.normal`
- :py:func:`~limix_plot.pca`
- :py:func:`~limix_plot.power`
- :py:func:`~limix_plot.qqplot`

.. _limix-plot: https://limix-plot.readthedocs.io/
"""
from __future__ import absolute_import

from limix_plot import (
    ConsensusCurve,
    box_aspect,
    get_pyplot,
    image,
    kinship,
    load_dataset,
    manhattan,
    normal,
    pca,
    power,
    qqplot,
    show,
)

__all__ = [
    "box_aspect",
    "ConsensusCurve",
    "image",
    "kinship",
    "load_dataset",
    "manhattan",
    "normal",
    "pca",
    "power",
    "qqplot",
    "image",
    "get_pyplot",
    "show",
]
