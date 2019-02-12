"""
io module
=========

Subpackage for I/O operations.

Modules
-------
bgen
csv
gen
hdf5
npy
plink
bimbam

Functions
---------
fetch
infer_filetype
"""

from . import bgen, bimbam, csv, gen, hdf5, npy, plink
from ._detect import infer_filetype
from ._fetch import fetch

__all__ = [
    "bgen",
    "csv",
    "gen",
    "hdf5",
    "npy",
    "plink",
    "bimbam",
    "fetch",
    "infer_filetype",
]
