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
from __future__ import absolute_import

__version__ = "2.0.0a3"

from . import display, glmm, her, io, plot, qc, qtl, sh, stats, example
from ._cli import cli
from ._config import config
from ._testit import test
from ._threads import set_max_nthreads, get_max_nthreads


__all__ = [
    "__version__",
    "test",
    "io",
    "plot",
    "qc",
    "qtl",
    "stats",
    "her",
    "glmm",
    "set_max_nthreads",
    "main",
    "sh",
    "config",
    "display",
    "cli",
    "get_max_nthreads",
    "example",
]
