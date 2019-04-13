"""
_data module
============

Private subpackage for standardisation of data structures and names.
"""
from ._asarray import asarray
from ._assert import assert_filetype, assert_target
from ._conf import CONF
from ._conform import conform_dataset
from ._lik import normalize_likelihood

__all__ = [
    "CONF",
    "asarray",
    "assert_filetype",
    "assert_target",
    "conform_dataset",
    "normalize_likelihood",
]
