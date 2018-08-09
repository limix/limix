from numpy import all as npall
from numpy import asarray, isfinite, ones
from pandas import DataFrame, Index
import dask.dataframe as da
import xarray as xr

from .likelihood import normalise_extreme_values


def _default_covariates_index(n):
    return ["covariate{}".format(i) for i in range(n)]


def _default_candidates_index(n):
    return ["candidate{}".format(i) for i in range(n)]


def _index_set_intersection(arrs):
    """Indices that are present in every indexed DataFrame/Series."""
    index_set = None
    for a in arrs:

        if (
            hasattr(a, "index")
            and isinstance(a.index, Index)
            or isinstance(a, xr.DataArray)
        ):

            i = set(asarray(_getindex(a)))

            if index_set is None:
                index_set = i
            else:
                index_set = index_set.intersection(i)

    if index_set is None:
        return set()
    return index_set


def _same_order_if_possible(index_set, arrs):
    for a in arrs:
        if index_set == set(a.index) and len(index_set) == len(a.index):
            return a.index.copy()
    return Index(index_set)


def _isindexed(a):
    return hasattr(a, "index") or isinstance(a, xr.DataArray)


def _getindex(a):
    if hasattr(a, "index"):
        return a.index
    if "samples" in a.dims:
        return a["samples"]
    return a["xsamples"]


def _infer_samples_index(arrs):
    if not isinstance(arrs, (list, tuple)):
        arrs = [arrs]

    ok = assert_compatible_samples_arrays(arrs)
    if not ok:
        msg = "The provided arrays are sample-wise incompatible."
        msg += " Please, check the number of rows."
        raise ValueError(msg)

    if len(arrs) == 0:
        return Index()

    iarrs = [a for a in arrs if _isindexed(a)]
    if len(iarrs) > 0:
        index_set = _index_set_intersection(arrs)
    else:
        index_set = _default_candidates_index(arrs[0].shape[0])

    return _same_order_if_possible(index_set, iarrs)


def assert_compatible_samples_arrays(arrs):
    """Check index and number of rows.

    The arrays are compatible if they have the same number of rows or
    if every array is indexed. If they have the same number of rows but some
    arrays are indexed and some not, the indexed arrays need to show the same
    index order as to be compatible to the non-indexed ones.
    """
    if not isinstance(arrs, (list, tuple)):
        arrs = [arrs]

    s = set([len(a) for a in arrs])
    if len(s) == 0:
        return True

    iarrs = [a for a in arrs if _isindexed(a)]
    if len(arrs) == len(iarrs):
        return True

    if len(s) == 1:
        if len(iarrs) == 0:
            return True
        # Make sure index have the same order.
        index = iarrs[0].index
        return all([all(index == a.index) for a in iarrs])

    return False


def _make_sure_phenotype_dataframe(y, samples_index):
    if _isindexed(y):
        y = y.loc[asarray(samples_index)]
    else:
        y = DataFrame(data=y, index=samples_index.copy())

    return y


def _make_sure_candidates_dataframe(G, samples_index):
    if _isindexed(G):
        G = G.loc[asarray(samples_index)]
    else:
        i = samples_index.copy()
        colnames = _default_candidates_index(G.shape[1])
        G = DataFrame(data=G, index=i, columns=colnames)

    return G


def _make_sure_covariates_dataframe(M, samples_index):
    if _isindexed(M):
        M = M.loc[asarray(samples_index)]
    else:
        i = samples_index.copy()
        colnames = _default_covariates_index(M.shape[1])
        M = DataFrame(data=M, index=i, columns=colnames)

    return M


def _make_sure_kinship_dataframe(K, samples_index):
    if _isindexed(K):
        K = K.loc[asarray(samples_index), :]
        K = K.loc[:, asarray(samples_index)]
    else:
        rows = samples_index.copy()
        cols = samples_index.copy()
        K = DataFrame(data=K, index=rows, columns=cols)

    return K


def _intersect_indices(arrs):
    indices = arrs[0].index.copy()
    for a in arrs:
        indices = indices.intersection(a.index)
    return indices


def _check_duplicity(samples_index, arrs):
    for a in arrs:
        indices = samples_index.intersection(asarray(_getindex(a)))
        if len(indices) > len(set(indices)):
            return False
    return True


def normalise_dataset(y, lik, M=None, G=None, K=None):
    """Convert to dataframes and align the arrays."""
    if isinstance(y, (tuple, list)):
        y = asarray(y, float).T

    if M is None:
        M = ones((y.shape[0], 1))
        if _isindexed(y):
            M = DataFrame(M, index=y.index)

    arrs = [a for a in [y, G, M, K] if a is not None]

    # Exit early if all the indices are the same.
    if all([_isindexed(a) for a in arrs]):
        a0 = arrs[0]
        if all([len(_getindex(a)) == len(_getindex(a0)) for a in arrs[1:]]):
            if all([all(_getindex(a) == _getindex(a0)) for a in arrs[1:]]):
                return dict(y=y, M=M, G=G, K=K)

    samples_index = _infer_samples_index(arrs)
    arrs0 = [a for a in [G, M, K] if a is not None and _isindexed(a)]
    ok = _check_duplicity(samples_index, arrs0)

    if not ok:
        msg = "Duplicated indices are allowed only on the"
        msg += " outcome variable when the indices differ."
        raise ValueError(msg)

    if len(samples_index) == 0:
        raise ValueError(
            "Could not infer an index for samples."
            " Please, check if the passed arrays are in order."
        )

    y = _make_sure_phenotype_dataframe(y, samples_index)
    samples_index = y.index

    if G is not None:
        G = _make_sure_candidates_dataframe(G, samples_index)
        if not isinstance(G, (da.DataFrame, xr.DataArray)) and not npall(isfinite(G)):
            raise ValueError("Candidate values must be finite.")

    M = _make_sure_covariates_dataframe(M, samples_index)

    if K is not None:
        K = _make_sure_kinship_dataframe(K, samples_index)
        if not isinstance(K, da.DataFrame) and not npall(isfinite(K)):
            raise ValueError("Covariance values must be finite.")

    y = normalise_extreme_values(y, lik)

    return dict(y=y, M=M, G=G, K=K)
