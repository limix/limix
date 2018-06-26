from __future__ import division

from numpy import clip, eye, ones, ascontiguousarray, atleast_2d
from numpy import stack as npy_stack
from numpy import concatenate, minimum
from pandas import DataFrame, Series, to_numeric, Int64Index, Index
from dask.array import Array as DaskArray
from dask.array import stack as dsk_stack
from dask.dataframe import DataFrame as DaskDataFrame

from .fprint import wprint
from .qc import normalise_covariance
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
    ncmax = 2 + (lik == "binomial")
    if len(y.columns) > ncmax:
        raise ValueError("The phenotype data frame has too many columns.")

    nc = 2 + (lik == "binomial")

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
        M = DataFrame({"offset": ones(nsamples, float)})
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

    if lik == "poisson":
        y = clip(y, 0., 25000.)

    y = y.astype(float)

    if not _all(isfinite(y)):
        msg = "One or more values of the provided phenotype "
        msg += "is not finite."
        raise ValueError(msg)

    if lik == "binomial":
        if y.ndim != 2:
            raise ValueError("Binomial phenotype matrix has to be" " bi-dimensional")
        if y.shape[1] != 2:
            raise ValueError("Binomial phenotype matrix has to have" " two columns.")

        v = y.values
        ratio = v[:, 0] / v[:, 1]
        v[:, 1] = minimum(v[:, 1], 300)
        v[:, 0] = ratio * v[:, 1]
        v[:, 0] = v[:, 0].round()
        y[:] = v
    return y


def kinship_process(K, nsamples, verbose):

    if K is None:
        K = eye(nsamples)

    K = asarray(K)

    ah = array_hash(K)
    nvalid = _cache["K"]["K"] is None
    nvalid = nvalid or (_cache["K"]["hash"] != ah)

    if nvalid:

        if not _all(isfinite(K)):
            msg = "One or more values of "
            msg += "the provided covariance matrix "
            msg += "is not finite."
            raise ValueError(msg)

        K = normalise_covariance(K)

        desc = "Eigen decomposition of the covariance matrix..."
        with Timer(desc=desc, disable=not verbose):
            QS = economic_qs(K)
            _cache["K"]["hash"] = ah
            _cache["K"]["QS"] = QS
            _cache["K"]["K"] = K
    else:
        QS = _cache["K"]["QS"]
        K = _cache["K"]["K"]

    return K, QS


def _isdataframe(df):
    return isinstance(df, (DataFrame, DaskDataFrame))


def _normalise_dataframe_phenotype(y, lik):
    if lik == "binomial":
        if y.shape[1] == 3:
            print(
                "Phenotype dataframe has three columns. I will be using the"
                " first one as the index."
            )
            y = y.set_index(y.columns[0])
    else:
        if y.shape[1] == 2:
            print(
                "Phenotype dataframe has two columns. I will be using the"
                " first one as the index."
            )
            y = y.set_index(y.columns[0])
    y = y.astype(float)
    return y


def normalise_phenotype_matrix(y, lik):
    if _isdataframe(y):
        return _normalise_dataframe_phenotype(y, lik)

    y = ascontiguousarray(y, float)
    y = atleast_2d(y.T).T

    if y.shape[1] > 2:
        raise ValueError("Outcome matrix must have two or one columns.")

    if lik == "binomial":
        pass
    else:
        if y.shape[1] > 1:
            y = DataFrame(data=y)
        else:
            index = default_samples_index(y.shape[0])
            y = DataFrame(data=y, index=index)

    return _normalise_dataframe_phenotype(y, lik)


def default_samples_index(n):
    return ["sample{}".format(i) for i in range(n)]


def default_covariates_index(n):
    return ["covariate{}".format(i) for i in range(n)]


def default_candidates_index(n):
    return ["candidate{}".format(i) for i in range(n)]


def _normalise_covariates_dataframe(M):
    M = M.astype(float)
    return M


def _normalise_candidates_dataframe(M):
    M = M.astype(float)
    return M


def normalise_covariates_matrix(M):
    if _isdataframe(M):
        return _normalise_covariates_dataframe(M)

    M = ascontiguousarray(M, float)
    M = atleast_2d(M.T).T

    index = default_samples_index(M.shape[0])
    columns = default_covariates_index(M.shape[1])
    M = DataFrame(M, columns=columns, index=index)

    return _normalise_covariates_dataframe(M)


def normalise_candidates_matrix(X):
    if _isdataframe(X):
        return _normalise_candidates_dataframe(X)

    X = ascontiguousarray(X, float)
    X = atleast_2d(X.T).T

    index = default_samples_index(X.shape[0])
    columns = default_candidates_index(X.shape[1])
    X = DataFrame(X, columns=columns, index=index)

    return _normalise_candidates_dataframe(X)


def _infer_same_index_col_type(M):
    a = M.index
    b = M.columns
    X = npy_stack((a.ravel(), b.ravel()), axis=0)
    X = Series(data=concatenate([a.ravel(), b.ravel()]))
    try:
        X = to_numeric(X).values
    except ValueError:
        X = X.astype("unicode").values.astype("unicode")

    dtype = X.dtype

    M.index = asarray(M.index, dtype=dtype)
    M.columns = asarray(M.columns, dtype=dtype)

    return M


def _normalise_kinship_dataframe(M):
    if M.shape[1] == M.shape[0] + 1:
        print(
            "Kinship matrix has one column too many. I will use the first"
            " one as the sample index."
        )
        M = M.set_index(M.columns[0])
        if isinstance(M.columns, Int64Index):
            M.columns = M.columns - 1

    elif M.shape[0] == M.shape[1] + 1:
        print(
            "Kinship matrix has one row too many. I will use the first"
            " one as the sample index."
        )
        M.columns = M.iloc[0]
        M = M.reindex(M.index.drop(0))
        if isinstance(M.index, Int64Index):
            M.index = M.index - 1

    M = M.astype(float)
    M = _infer_same_index_col_type(M)

    return M


def normalise_kinship_matrix(M):
    if _isdataframe(M):
        return _normalise_kinship_dataframe(M)

    M = ascontiguousarray(M, float)
    M = atleast_2d(M.T).T

    index = default_samples_index(M.shape[0])
    columns = default_samples_index(M.shape[1])
    M = DataFrame(M, columns=columns, index=index)

    return _normalise_kinship_dataframe(M)


def _index_set_intersection(arrs):
    index_set = None
    for a in arrs:

        if hasattr(a, "index") and isinstance(a.index, Index):

            i = set(a.index)
            if len(i) < len(a.index):
                raise ValueError("Duplicated index.")

            if index_set is None:
                index_set = i
            else:
                index_set = index_set.intersection(i)

    if index_set is None:
        return set()
    return index_set


def _equal_index_if_possible(index_set, arrs):
    for a in arrs:
        if hasattr(a, "index"):
            if index_set == set(a.index):
                return a.index.copy()
    return Index(index_set)


def infer_samples_index(arrs):
    index_set = _index_set_intersection(arrs)
    if len(index_set) == 0:
        index_set = default_candidates_index(arrs[0].shape[0])
    return _equal_index_if_possible(index_set, arrs)
