from __future__ import absolute_import

from . import bgen, bimbam, csv, gen, hdf5, npy, plink
from ._phenotype import fetch_phenotype
from ._genotype import fetch_genotype
from ._detect import detect_file_type

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
    "detect_file_type",
]
