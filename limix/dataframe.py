from numpy import all as npall
from numpy import asarray, isfinite, ones
from pandas import DataFrame, Index

from .likelihood import normalise_extreme_values


def _default_covariates_index(n):
    return ['covariate{}'.format(i) for i in range(n)]


def _default_candidates_index(n):
    return ['candidate{}'.format(i) for i in range(n)]


def _index_set_intersection(arrs):
    """Return the indices that are present in every indexed DataFrame/Series.

    It fails if there are duplicated indices.
    """
    index_set = None
    for a in arrs:

        if hasattr(a, 'index') and isinstance(a.index, Index):

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


def _same_order_if_possible(index_set, arrs):
    for a in arrs:
        if index_set == set(a.index):
            return a.index.copy()
    return Index(index_set)


def _isindexed(a):
    return hasattr(a, 'index') and isinstance(a.index, Index)


def infer_samples_index(arrs):
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
        y = y.loc[y.index.intersection(samples_index)]
    else:
        y = DataFrame(data=y, index=samples_index.copy())

    return y


def _make_sure_candidates_dataframe(G, samples_index):
    if _isindexed(G):
        G = G.loc[G.index.intersection(samples_index)]
    else:
        i = samples_index.copy()
        colnames = _default_candidates_index(G.shape[1])
        G = DataFrame(data=G, index=i, columns=colnames)

    return G


def _make_sure_covariates_dataframe(M, samples_index):
    if _isindexed(M):
        M = M.loc[M.index.intersection(samples_index)]
    else:
        i = samples_index.copy()
        colnames = _default_covariates_index(M.shape[1])
        M = DataFrame(data=M, index=i, columns=colnames)

    return M


def _make_sure_kinship_dataframe(K, samples_index):
    if _isindexed(K):
        K = K.loc[K.index.intersection(samples_index), :]
        K = K.loc[:, K.columns.intersection(samples_index)]
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


def normalise_dataset(y, lik, M=None, G=None, K=None):
    y = asarray(y, float).T

    if M is None:
        M = ones((y.shape[0], 1))

    arrs = [a for a in [y, G, M, K] if a is not None]
    samples_index = infer_samples_index(arrs)
    if len(samples_index) == 0:
        raise ValueError("Could not infer an index for samples."
                         " Please, check if the passed arrays are in order.")

    y = _make_sure_phenotype_dataframe(y, samples_index)
    if G is not None:
        G = _make_sure_candidates_dataframe(G, samples_index)
        if not npall(isfinite(G)):
            raise ValueError("Candidate values must be finite.")

    M = _make_sure_covariates_dataframe(M, samples_index)

    if K is not None:
        K = _make_sure_kinship_dataframe(K, samples_index)

    y = normalise_extreme_values(y, lik)

    indices = _intersect_indices([a for a in [y, G, M, K] if a is not None])

    M = M.loc[indices]
    if G is not None:
        G = G.loc[indices]
    y = y.loc[indices]

    if K is not None:
        K = K.loc[indices, :].loc[:, indices]
        if not npall(isfinite(K)):
            raise ValueError("Kinship matrix contain non-finite values.")

    return dict(y=y, M=M, G=G, K=K)
