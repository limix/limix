from __future__ import absolute_import

from . import bgen, bimbam, csv, hdf5, npy, plink
from . import gen
from ._fetch import fetch
from ._detect import infer_filetype

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
