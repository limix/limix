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
from .detect import detect_file_type

__all__ = ["bgen", "csv", "examples", "detect_file_type", "gen", "hdf5", "npy", "plink"]
