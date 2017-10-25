from __future__ import division

from numpy import clip, eye, ones
from numpy import stack as npy_stack
from dask.array import Array as DaskArray
from dask.array import stack as dsk_stack
from dask.dataframe import DataFrame as DaskDataFrame

from .qc import gower_norm
from .util import Timer, array_hash
from .util.npy_dask import all as _all
from .util.npy_dask import asarray, isfinite
from numpy_sugar.linalg import economic_qs

_cache = dict(K=dict(K=None, QS=None, hash=None))


def assure_named_matrix(A):
    from pandas import DataFrame

    if isinstance(A, DataFrame):
        return A

    if isinstance(A, DaskDataFrame):
        return A

    A = asarray(A)
    if A.ndim != 2:
        msg = "Wrong number of dimensions. "
        msg += "It should be two."
        raise ValueError(msg)
    A = DataFrame(A, columns=list(range(A.shape[1])))

    return A


def named_to_unamed_matrix(M):
    from pandas import Series

    k = next(iter(M.keys()))

    if isinstance(M[k], Series):
        return npy_stack(list(M.values), axis=0)

    if isinstance(M[k], DaskArray):
        return dsk_stack(list(M.values), axis=1)

    return npy_stack(list(M.values), axis=1)


def covariates_process(M, nsamples):
    from pandas import DataFrame

    if M is None:
        M = DataFrame({'offset': ones(nsamples, float)})
    else:
        M = assure_named_matrix(M)

    if not _all(isfinite(M.values)):
        msg = "One or more values of the provided covariates "
        msg += "is not finite."
        raise ValueError(msg)

    if M.shape[0] != nsamples:
        raise ValueError("Wrong number of rows in the covariate matrix.")

    return M


def phenotype_process(lik, y):
    if lik == 'poisson':
        y = clip(y, 0., 25000.)

    if isinstance(y, (tuple, list)):
        y = tuple([asarray(p, float) for p in y])
    else:
        y = asarray(y, float)

    if not _all(isfinite(y)):
        msg = "One or more values of the provided phenotype "
        msg += "is not finite."
        raise ValueError(msg)
    return y


def kinship_process(K, nsamples, verbose):

    if K is None:
        K = eye(nsamples)

    K = asarray(K)

    ah = array_hash(K)
    nvalid = _cache['K']['K'] is None
    nvalid = nvalid or (_cache['K']['hash'] != ah)

    if nvalid:

        if not _all(isfinite(K)):
            msg = "One or more values of "
            msg += "the provided covariance matrix "
            msg += "is not finite."
            raise ValueError(msg)

        K = gower_norm(K)

        desc = "Eigen decomposition of the covariance matrix..."
        with Timer(desc=desc, disable=not verbose):
            QS = economic_qs(K)
            _cache['K']['hash'] = ah
            _cache['K']['QS'] = QS
            _cache['K']['K'] = K
    else:
        QS = _cache['K']['QS']
        K = _cache['K']['K']

    return K, QS
