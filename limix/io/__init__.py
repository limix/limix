r"""
**********
I/O module
**********

PLINK reader
^^^^^^^^^^^^

.. autofunction:: limix.io.read_plink

HDF5 reader
^^^^^^^^^^^

.. :autoclass:: limix.io.h5data_fetcher

CSV reader
^^^^^^^^^^

.. autofunction:: limix.io.read_csv

GEN reader
^^^^^^^^^^

.. autofunction:: limix.io.read_gen

"""

from . import examples
from ._csv import read_csv
from .gen import read_gen
from . import hdf5
from .hdf5 import h5data_fetcher
from .plink import read_plink
from .util import file_type
from . import _csv as csv
from . import plink
from . import npy


def genotype_reader(*args, **kwargs):
    from limix_legacy.io import genotype_reader
    return genotype_reader(*args, **kwargs)


def phenotype_reader(*args, **kwargs):
    from limix_legacy.io import phenotype_reader
    return phenotype_reader(*args, **kwargs)


def data(*args, **kwargs):
    from limix_legacy.io import data
    return data(*args, **kwargs)


__all__ = [
    'read_plink', 'h5data_fetcher', 'read_csv', 'read_gen', 'examples',
    'genotype_reader', 'phenotype_reader', 'data', 'file_type', 'hdf5',
    'csv', 'plink', 'npy'
]
