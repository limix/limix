from __future__ import absolute_import

from . import bgen, bimbam, csv, hdf5, npy, plink
from . import gen
from ._phenotype import fetch_phenotype
from ._genotype import fetch_genotype
from ._detect import detect_filetype, get_fetch_spec

__all__ = [
    "bgen",
    "csv",
    "gen",
    "hdf5",
    "npy",
    "plink",
    "bimbam",
    "fetch_phenotype",
    "fetch_genotype",
    "detect_filetype",
    "get_fetch_spec",
]
