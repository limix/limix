r"""
*************
limix package
*************

A flexible and fast mixed model toolbox.

"""

from __future__ import absolute_import as _

from . import cmdlimix, heritability, io, qc, qtl, stats, util
from .testit import test
from .threads import get_max_nthreads, set_max_nthreads
from .util import download, extract, filehash
from .cmd import call
from ._examples import load_dataset
from ._plot import plot, ConsensusCurve, LimixPlot

__name__ = "limix"
__version__ = "1.1.0"
__author__ = "C. Lippert, D. Horta, F. P. Casale, and O. Stegle"
__author_email__ = "stegle@ebi.ac.uk"

__all__ = [
    "__name__", "__version__", "__author__", "__author_email__", "test", 'io',
    'plot', 'qc', 'qtl', 'stats', 'util', 'heritability', 'download',
    'extract', 'filehash', 'cmdlimix', 'set_max_nthreads', 'get_max_nthreads',
    'call', 'ConsensusCurve', 'LimixPlot', 'load_dataset'
]
