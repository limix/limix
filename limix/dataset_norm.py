from numpy import unique, asarray
from xarray import DataArray


def normalise_dataset(y, M, G=None, K=None):

    y = DataArray(y, encoding={"dtype": "float64"})
    y = y.rename({y.dims[0]: "sample"})

    idx = y.coords["sample"].values
    if len(unique(idx)) < len(idx):
        raise ValueError("Non-unique sample ids are not allowed in the outcome.")

    M = DataArray(M, encoding={"dtype": "float64"})
    M = M.rename({M.dims[0]: "sample"})

    if G is not None:
        G = DataArray(G, encoding={"dtype": "float64"})
        G = G.rename({G.dims[0]: "sample"})

    if K is not None:
        K = DataArray(K, encoding={"dtype": "float64"})
        K = K.rename({K.dims[0]: "sample_0", K.dims[1]: "sample_1"})

    arrs = [a for a in [y, M, G] if a is not None]
    if len(y.dims) == 1:
        y = y.expand_dims("trait", 1)

    if K is not None:
        arrs += [K, K.T]

    M = M.sel(sample=y.coords["sample"].values)
    K = K.sel(sample_0=y.coords["sample"].values)
    K = K.sel(sample_1=y.coords["sample"].values)
    if G is not None:
        G = G.sel(sample=y.coords["sample"].values)

    y = y.rename({y.dims[1]: "trait"})
    M = M.rename({M.dims[1]: "covariate"})
    if G is not None:
        G = G.rename({G.dims[1]: "candidate"})

    return dict(y=y, M=M, G=G, K=K)


def _create_index(vals):
    a = asarray(vals, object)
    b = asarray(vals, object)
    return DataArray(a, dims=["sample"], coords={"sample": b})


def _infer_samples_index(arrs):
    ok = _check_samples_arrays_compatibility(arrs)
    if not ok:
        msg = "The provided arrays are sample-wise incompatible."
        msg += " Please, check the number of rows."
        raise ValueError(msg)

    if len(arrs) == 0:
        return _create_index([])

    iarrs = [a for a in arrs if _has_sample_index(a)]
    if len(iarrs) == 0:
        return _create_index(range(arrs[0].shape[0]))

    index_set = _index_set_intersection(iarrs)

    return _create_index(_same_order_if_possible(index_set, iarrs))


def _check_samples_arrays_compatibility(arrs):
    """Check index and number of rows.

    The arrays are compatible if they have the same number of rows or
    if every array is indexed. If they have the same number of rows but some
    arrays are indexed and some not, the indexed arrays need to show the same
    index order as to be compatible to the non-indexed ones.
    """
    s = set([len(a) for a in arrs])
    if len(s) == 0:
        return True

    iarrs = [a for a in arrs if _has_sample_index(a)]
    if len(arrs) == len(iarrs):
        return True

    if len(s) == 1:
        if len(iarrs) == 0:
            return True
        # Make sure index have the same order.
        index = _get_sample_index(iarrs[0])
        return all([all(index == a.index) for a in iarrs])

    return False


def _has_sample_index(x):
    if hasattr(x, "index"):
        return True
    return hasattr(x, "coords") and "sample" in x.dims


def _get_sample_index(x):
    if hasattr(x, "index"):
        return x.index.values
    return x.coords["sample"].values


def _index_set_intersection(arrs):
    """Indices that are present in every indexed DataFrame/Series."""
    index_set = None
    for a in arrs:
        i = set(asarray(_get_sample_index(a)))

        if index_set is None:
            index_set = i
        else:
            index_set = index_set.intersection(i)

    if index_set is None:
        return set()

    return index_set


def _same_order_if_possible(index_set, arrs):
    for a in arrs:
        i = _get_sample_index(a)
        if len(index_set) == len(i) and index_set == set(i):
            return i.copy()
    return _create_index([])
