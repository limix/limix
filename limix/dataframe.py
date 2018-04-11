from pandas import DataFrame, Index


def make_sure_dataframes(y, G=None, M=None, K=None):

    # Rows represent samples in every array of this list.
    arrs = [a for a in [y, G, M, K] if a is not None]

    samples_index = infer_samples_index(arrs)
    # if len(samples_index) == 0:
    #     raise ValueError("Could not infer an index for samples."
    #                      " Please, check if the passed arrays are in order.")
    #
    # if hasattr(y, 'index'):
    #     y = y.loc[y.index.intersection(samples_index)]
    # else:
    #     y = DataFrame(data=y, index=samples_index.copy())
    #
    # if hasattr(G, 'index'):
    #     G = G.loc[G.index.intersection(samples_index)]
    # else:
    #     G = DataFrame(
    #         data=G,
    #         index=samples_index.copy(),
    #         columns=default_candidates_index(G.shape[1]))
    #
    # if hasattr(M, 'index'):
    #     M = M.loc[M.index.intersection(samples_index)]
    # else:
    #     M = DataFrame(
    #         data=M,
    #         index=samples_index.copy(),
    #         columns=default_covariates_index(M.shape[1]))
    #
    # if K is not None:
    #     if hasattr(K, 'index'):
    #         K = K.loc[K.index.intersection(samples_index), :]
    #         K = K.loc[:, K.columns.intersection(samples_index)]
    #     else:
    #         K = DataFrame(
    #             data=K,
    #             index=samples_index.copy(),
    #             columns=samples_index.copy())
    #
    # y = normalise_phenotype_matrix(y, lik)
    # if K is not None:
    #     K = normalise_kinship_matrix(K)
    #
    # G = normalise_candidates_matrix(G)
    # M = covariates_process(M, nsamples)
    # M = normalise_covariates_matrix(M)
    #
    # y = phenotype_process(lik, y)
    #
    # if not npall(isfinite(G)):
    #     raise ValueError("Variant values must be finite.")
    #
    # mixed = K is not None
    #
    # indices = _intersect_indices(G, y, K, M)
    #
    # G = G.loc[indices, :]
    # y = y.loc[indices, :]
    #
    # if K is not None:
    #     K = K.loc[indices, :]
    #     K = K.loc[:, indices]
    #
    # M = M.loc[indices, :]
    return dict(y=y, G=G, M=M, K=K)


def default_samples_index(n):
    return ['sample{}'.format(i) for i in range(n)]


def default_covariates_index(n):
    return ['covariate{}'.format(i) for i in range(n)]


def default_candidates_index(n):
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


def _equal_index_if_possible(index_set, arrs):
    for a in arrs:
        if hasattr(a, 'index'):
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
        index_set = default_candidates_index(arrs[0].shape[0])

    # if len(index_set) == 0:
    #     index_set = default_candidates_index(arrs[0].shape[0])
    return _equal_index_if_possible(index_set, arrs)


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
