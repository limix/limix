from __future__ import division

from numpy import clip, eye, ones
from numpy import stack as npy_stack
from pandas import DataFrame, Series
from dask.array import Array as DaskArray
from dask.array import stack as dsk_stack
from dask.dataframe import DataFrame as DaskDataFrame

from .fprint import wprint
from .qc import gower_norm
from .util import Timer, array_hash
from .util.npy_dask import all as _all
from .util.npy_dask import asarray, isfinite
from numpy_sugar.linalg import economic_qs

_cache = dict(K=dict(K=None, QS=None, hash=None))


def assure_named_phenotype(lik, y):
    if isinstance(y, DataFrame):
        return _normalise_phenotype_dataframe(lik, y)

    if isinstance(y, DaskDataFrame):
        return _normalise_phenotype_dataframe(lik, y)

    y = asarray(y, float)
    y = DataFrame(y)
    y = _normalise_phenotype_dataframe(lik, y)

    return y


def _normalise_phenotype_dataframe(lik, y):
    ncmax = 2 + (lik == 'binomial')
    if len(y.columns) > ncmax:
        raise ValueError("The phenotype data frame has too many columns.")

    nc = 2 + (lik == 'binomial')

    if len(y.columns) == nc:
        msg = "The ({}) phenotype data frame has {} columns. ".format(lik, nc)
        msg += "I will consider that"
        msg += " the first one is a column of row indices."
        wprint(msg)
        y = y.set_index(y.columns[0])
        y.index.name = None

    if len(y.columns) == 0:
        raise ValueError("The phenotype data frame has no columns.")

    return y


def assure_named_matrix(A):
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


def assure_named_columns(A):
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
    k = next(iter(M.keys()))

    if isinstance(M[k], Series):
        return npy_stack(list(M.values), axis=0)

    if isinstance(M[k], DaskArray):
        return dsk_stack(list(M.values), axis=1)

    return npy_stack(list(M.values), axis=1)


def covariates_process(M, nsamples):
    if M is None:
        M = DataFrame({'offset': ones(nsamples, float)})
    else:
        M = assure_named_columns(M)

    if not _all(isfinite(M.values)):
        msg = "One or more values of the provided covariates "
        msg += "is not finite."
        raise ValueError(msg)

    if M.shape[0] != nsamples:
        raise ValueError("Wrong number of rows in the covariate matrix.")

    return M


def phenotype_process(lik, y):
    if isinstance(y, (tuple, list)):
        y = asarray(y, float).T
    y = assure_named_phenotype(lik, y)

    if lik == 'poisson':
        y = clip(y, 0., 25000.)

    y = y.astype(float)

    if not _all(isfinite(y)):
        msg = "One or more values of the provided phenotype "
        msg += "is not finite."
        raise ValueError(msg)

    if lik == 'binomial':
        if y.ndim != 2:
            raise ValueError("Binomial phenotype matrix has to be"
                             " bi-dimensional")
        if y.shape[1] != 2:
            raise ValueError("Binomial phenotype matrix has to have"
                             " two columns.")
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
