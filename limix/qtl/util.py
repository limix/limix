from __future__ import division

from collections import OrderedDict

import dask.array as da
from numpy import stack as npy_stack
from numpy import ones


def assure_named(M, nsamples):
    from pandas import DataFrame
    from limix.util.npy_dask import asarray

    if isinstance(M, DataFrame):
        return M

    M = asarray(M)
    M = DataFrame(M, columns=list(range(M.shape[1])))

    return M
