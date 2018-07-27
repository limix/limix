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

from . import cmdlimix, her, io, qc, qtl, stats, util, plot, glmm
from .cmd import call
from .cmdlimix import main
from .limits import is_large_file
from ._testit import test
from .threads import get_max_nthreads, set_max_nthreads
from .util import download, extract, filehash

__version__ = "2.0.0-alpha.1"

__all__ = [
    "__version__",
    "test",
    "io",
    "plot",
    "qc",
    "qtl",
    "stats",
    "util",
    "her",
    "glmm",
    "download",
    "extract",
    "filehash",
    "cmdlimix",
    "set_max_nthreads",
    "get_max_nthreads",
    "call",
    "is_large_file",
    "normalise_phenotype",
    "main",
]
