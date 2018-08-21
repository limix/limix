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
from __future__ import absolute_import

from . import bgen, examples, hdf5, npy, plink, gen, csv
from .detect import file_type, possible_file_types

__all__ = [
    "examples",
    "file_type",
    "hdf5",
    "csv",
    "plink",
    "npy",
    "possible_file_types",
    "bgen",
    "gen",
]
