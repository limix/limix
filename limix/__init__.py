r"""
*************
limix package
*************

A flexible and fast mixed model toolbox.

"""

from __future__ import absolute_import as _

from . import (
    heritability, io, iset, mtset, plot, qc, qtl, scripts, stats, util, vardec
)

__name__ = "limix"
__version__ = "1.1.0"
__author__ = ("Christoph Lippert, Danilo Horta," +
              " Francesco Paolo Casale, and Oliver Stegle")
__author_email__ = "stegle@ebi.ac.uk"

__all__ = [
    "__name__", "__version__", "__author__", "__author_email__", "test", 'io',
    'plot', 'qc', 'qtl', 'stats', 'util', 'vardec', 'mtset', 'iset', 'scripts',
    'heritability'
]
