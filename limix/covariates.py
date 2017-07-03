from __future__ import division

from collections import OrderedDict

import dask.array as da
from numpy import stack as npy_stack
from numpy import ones

from limix.util import asarray


def assure_named_covariates(M, nsamples):

    if isinstance(M, dict):
        return M

    if M is None:
        M = OrderedDict(offset=ones(nsamples, float))
    else:
        M = asarray(M)
        M = [('C%d' % i, M[:, i]) for i in range(M.shape[1])]
        M = OrderedDict(M)

    return M


def named_covariates_to_array(M):
    k = next(iter(M.keys()))

    if isinstance(M[k], da.Array):
        return da.stack(list(M.values()), axis=1)

    return npy_stack(list(M.values()), axis=1)
