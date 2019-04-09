"""
_bits module
============

Private subpackage for low-level operations.

Modules
-------
dask
deco
pandas
xarray

Functions
---------
get_shape
"""
from . import dask, deco, numpy, pandas, xarray
from .array import cdot, get_shape, unvec, vec

__all__ = [
    "dask",
    "xarray",
    "deco",
    "pandas",
    "get_shape",
    "numpy",
    "vec",
    "unvec",
    "cdot",
]
