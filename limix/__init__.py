"""
Limix package
=============

A flexible and fast generalised mixed model toolbox.

The official documentation together with examples and tutorials can be found
at https://limix.readthedocs.io/.
"""
from . import __config__, her, io, plot, qc, qtl, sh, stats, threads, vardec
from ._cli import cli
from ._config import config
from ._example import file_example
from ._testit import test

__version__ = "3.0.2"


__all__ = [
    "__config__",
    "__version__",
    "cli",
    "config",
    "file_example",
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
    "vardec",
    "model",
]
