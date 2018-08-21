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

from . import bgen, csv, examples, gen, hdf5, npy, plink
from .detect import file_type, possible_file_types

__all__ = [
    "bgen",
    "csv",
    "examples",
    "file_type",
    "gen",
    "hdf5",
    "npy",
    "plink",
    "possible_file_types",
]
