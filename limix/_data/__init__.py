"""
_data package
=============

Private subpackage for standardisation of data structures and names.

Modules
-------
_conform    Standardise data structures.
_data       Manipulation of phenotype, genotype, and other data names.
_dataarray  Manipulation of names stored in DataArray.
_dim        DataArray dimension name manipulation.

"""
from __future__ import absolute_import as _

from ._conform import conform_dataset, to_dataarray
from ._dataarray import rename_dims
from ._dim import (
    dim_hint_to_name,
    dim_name_to_hint,
    get_dim_axis,
    get_dims_from_data_name,
    is_dim_hint,
    is_dim_name,
)
from ._data import (
    get_short_data_names,
    is_data_name,
    is_short_data_name,
    to_data_name,
    to_short_data_name,
)
from ._lik import check_likelihood_name

__all__ = [
    "check_likelihood_name",
    "conform_dataset",
    "dim_hint_to_name",
    "dim_name_to_hint",
    "get_dim_axis",
    "get_dims_from_data_name",
    "get_short_data_names",
    "is_data_name",
    "is_dim_hint",
    "is_dim_name",
    "is_short_data_name",
    "rename_dims",
    "to_data_name",
    "to_dataarray",
    "to_short_data_name",
]
