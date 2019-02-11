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
from ._asarray import asarray
from ._conform import conform_dataset, to_dataarray
from ._dataarray import rename_dims
from ._dim import is_dim_hint
from ._lik import check_likelihood_name

__all__ = [
    "asarray",
    "conform_dataset",
    "to_dataarray",
    "rename_dims",
    "check_likelihood_name",
    "is_dim_hint",
]
