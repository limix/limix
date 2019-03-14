"""
_data module
============

Private subpackage for standardisation of data structures and names.
"""
from ._asarray import asarray
from ._assert import assert_filetype, assert_likelihood, assert_target
from ._conf import CONF
from ._conform import conform_dataset

__all__ = [
    "asarray",
    "conform_dataset",
    "CONF",
    "assert_target",
    "assert_filetype",
    "assert_likelihood",
]
