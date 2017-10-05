r"""
*************
limix package
*************

A flexible and fast mixed model toolbox.

"""

from __future__ import absolute_import as _

from . import heritability, io, plot, qc, qtl, scripts, stats, util
from .testit import test
from .util import download, extract, filehash

__name__ = "limix"
__version__ = "1.1.0"
__author__ = "C. Lippert, D. Horta, F. P. Casale, and O. Stegle"
__author_email__ = "stegle@ebi.ac.uk"

__all__ = [
    "__name__", "__version__", "__author__", "__author_email__", "test", 'io',
    'plot', 'qc', 'qtl', 'stats', 'util', 'scripts', 'heritability',
    'download', 'extract', 'filehash'
]
