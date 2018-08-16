from numpy import unique, asarray
from xarray import DataArray


def normalise_dataset(y, M, G=None, K=None):

    y = DataArray(y, encoding={"dtype": "float64"})
    y = y.rename({y.dims[0]: "sample"})

    M = DataArray(M, encoding={"dtype": "float64"})
    M = M.rename({M.dims[0]: "sample"})

    idx = M.coords["sample"].values
    if len(unique(idx)) < len(idx):
        msg = "Non-unique sample ids are not allowed in the covariates array."
        raise ValueError(msg)

    if G is not None:
        G = DataArray(G, encoding={"dtype": "float64"})
        G = G.rename({G.dims[0]: "sample"})

        idx = G.coords["sample"].values
        if len(unique(idx)) < len(idx):
            msg = "Non-unique sample ids are not allowed in the candidates array."
            raise ValueError(msg)

    if K is not None:
        K = DataArray(K, encoding={"dtype": "float64"})
        K = K.rename({K.dims[0]: "sample_0", K.dims[1]: "sample_1"})

        idx0 = K.coords["sample_0"].values
        idx1 = K.coords["sample_1"].values
        if len(unique(idx0)) < len(idx0) or len(unique(idx1)) < len(idx1):
            msg = "Non-unique sample ids are not allowed in the covariance matrix."
            raise ValueError(msg)

    arrs = [a for a in [y, M, G] if a is not None]
    if len(y.dims) == 1:
        y = y.expand_dims("trait", 1)

    if K is not None:
        arrs += [K, K.T]

    inc_msg = "The provided outcome and {} arrays are sample-wise incompatible."
    try:
        M = M.sel(sample=y.coords["sample"].values)
    except IndexError as e:
        msg = inc_msg.format("covariates")
        raise ValueError(str(e) + "\n\n" + inc_msg)

    if K is not None:
        try:
            K = K.sel(sample_0=y.coords["sample"].values)
            K = K.sel(sample_1=y.coords["sample"].values)
        except IndexError as e:
            msg = inc_msg.format("covariance")
            raise ValueError(str(e) + "\n\n" + inc_msg)

    if G is not None:
        try:
            G = G.sel(sample=y.coords["sample"].values)
        except IndexError as e:
            msg = inc_msg.format("candidates")
            raise ValueError(str(e) + "\n\n" + inc_msg)

    y = y.rename({y.dims[1]: "trait"})
    M = M.rename({M.dims[1]: "covariate"})
    if G is not None:
        G = G.rename({G.dims[1]: "candidate"})

    nsamples = len(y)
    if nsamples != M.shape[0]:
        msg = inc_msg.format("covariates")
        raise ValueError(inc_msg)

    if K is not None:
        if nsamples != K.shape[0] or nsamples != K.shape[1]:
            inc_msg = inc_msg.format("covariance")
            raise ValueError(inc_msg)

    if G is not None:
        if nsamples != G.shape[0]:
            inc_msg = inc_msg.format("candidates")
            raise ValueError(inc_msg)

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
