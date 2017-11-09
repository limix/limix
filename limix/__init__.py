r"""
limix package
=============

A flexible and fast generalised mixed model toolbox.

Modules
-------
cmd
    Command line interface.
her
    Genetic heritability estimation.
io
    Functions for reading common files used in genetics.
plot
    Visualization of data and results for genetic analysis.
qc
    Quality control for genetic data sets.
qtl
    Quantitative trait locus analysis.
stats
    PCA, confusion matrix, p-value correction, and others.

The official documentation together with examples and tutorials can be found
at https://limix.readthedocs.io/.
"""

from __future__ import absolute_import as _

from . import cmdlimix, her, io, qc, qtl, stats, util
from .testit import test
from .threads import get_max_nthreads, set_max_nthreads
from .util import download, extract, filehash
from .cmd import call
from ._examples import load_dataset
from ._plot import plot, ConsensusCurve, LimixPlot
from .limits import is_large_file

__name__ = "limix"
__version__ = "1.2.0"
__author__ = "C. Lippert, D. Horta, F. P. Casale, and O. Stegle"
__author_email__ = "stegle@ebi.ac.uk"

__all__ = [
    "__name__", "__version__", "__author__", "__author_email__", "test", 'io',
    'plot', 'qc', 'qtl', 'stats', 'util', 'her', 'download', 'extract',
    'filehash', 'cmdlimix', 'set_max_nthreads', 'get_max_nthreads', 'call',
    'ConsensusCurve', 'LimixPlot', 'load_dataset', 'is_large_file',
    'normalise_phenotype'
]
