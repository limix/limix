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

from . import _csv as csv
from . import bgen, examples, hdf5, npy, plink
from ._csv import read_csv
from .detect import file_type, possible_file_types
from .gen import read_gen
from .hdf5 import h5data_fetcher, read_hdf5_limix
from .plink import read_plink

__all__ = [
    'read_plink', 'h5data_fetcher', 'read_csv', 'read_gen', 'examples',
    'file_type', 'hdf5', 'csv', 'plink', 'npy', 'possible_file_types',
    'read_hdf5_limix', 'bgen'
]
