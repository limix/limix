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
gwas
    Methods for GWAS and QTL mapping.
stats
    PCA, confusion matrix, p-value correction, and others.

The official documentation together with examples and tutorials can be found
at https://limix.readthedocs.io/.
"""
from __future__ import absolute_import as _

from . import example, glmm, her, io, plot, qc, qtl, sh, stats, threads
from ._cli import cli
from ._config import config
from ._testit import test

__version__ = "2.0.6"


__all__ = [
    "__version__",
    "cli",
    "config",
    "example",
    "glmm",
    "gwas",
    "her",
    "io",
    "main",
    "plot",
    "qc",
    "qtl",
    "sh",
    "stats",
    "test",
    "threads",
]
