"""
_data module
============

Private subpackage for standardisation of data structures and names.

Functions
-------
conform_dataset       Standardise arrays jointly into limix DataArray.
asarray               Convert array-like to limix DataArray.
check_likelihood_name Check if a likelihood name is supported.
"""
from ._asarray import asarray
from ._assert import assert_filetype, assert_target
from ._conf import CONF
from ._conform import conform_dataset
from ._lik import check_likelihood_name

__all__ = [
    "asarray",
    "conform_dataset",
    "check_likelihood_name",
    "CONF",
    "assert_target",
    "assert_filetype",
]
