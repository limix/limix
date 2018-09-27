from collections import Counter
from numpy import unique, asarray
import dask.dataframe as dd
import dask.array as da
import xarray as xr
from xarray import DataArray
import numpy as np
from numpy import array
from ._dask import array_shape_reveal


def normalise_dataset(y, M, G=None, K=None):

    y = _dataarray_upcast(y)
    y = y.rename({y.dims[0]: "sample"})
    y = y.rename({y.dims[1]: "trait"})

    M = _dataarray_upcast(M)
    M = M.rename({M.dims[0]: "sample"})
    M = M.rename({M.dims[1]: "covariate"})

    if G is not None:
        G = _dataarray_upcast(G)
        G = G.rename({G.dims[0]: "sample"})
        G = G.rename({G.dims[1]: "candidate"})

    if K is not None:
        K = _dataarray_upcast(K)
        K = K.rename({K.dims[0]: "sample_0"})
        K = K.rename({K.dims[1]: "sample_1"})

    arrs = [a for a in [y, M, G, K] if a is not None]
    valid_samples, invalid_samples = _infer_samples_index(arrs)

    if "sample" not in y.coords:
        y = y.assign_coords(sample=list(valid_samples))
    else:
        # idx = y.get_index("sample")
        # dup_idx = set(idx[idx.duplicated()].unique())
        # uni_idx = set(idx) - dup_idx
        # R0 = idx.isin(dup_idx)
        # uni_idx & valid_samples
        # Index(uni_idx).sel(Counter(uni_idx & set(valid_samples)))
        # idx.isin(dup_idx & idx)
        # valid_samples
        y = y.loc[{"sample": y.get_index("sample").isin(valid_samples)}]

    if "sample" not in M.coords:
        M = M.assign_coords(sample=list(valid_samples))
    else:
        M = M.loc[{"sample": M.get_index("sample").isin(valid_samples)}]

    if G is not None:
        if "sample" not in G.coords:
            G = G.assign_coords(sample=list(valid_samples))
        else:
            G = G.loc[{"sample": G.get_index("sample").isin(valid_samples)}]

    if K is not None:
        for k in ["sample_0", "sample_1"]:
            if k not in K.coords:
                K = K.assign_coords(**{k: list(valid_samples)})
            else:
                K = K.loc[{k: K.get_index(k).isin(valid_samples)}]

    arrs = [a for a in [y, M, G] if a is not None]

    if K is not None:
        arrs += [K, K.T]

    idx = y.coords["sample"]

    for a in arrs:
        if len(a.coords[a.dims[0]]) != len(idx):
            break
    else:
        return dict(y=y, M=M, G=G, K=K)

    idx = M.coords["sample"].values
    if len(unique(idx)) < len(idx):
        msg = "Non-unique sample ids are not allowed in the covariates array"
        msg += " if the sample ids are not equal nor in the same order."
        raise ValueError(msg)

    if y.name is None:
        y.name = "outcome"

    if G is not None:
        idx = G.coords["sample"].values
        if len(unique(idx)) < len(idx):
            msg = "Non-unique sample ids are not allowed in the candidates array"
            msg += " if the sample ids are not equal nor in the same order."
            raise ValueError(msg)
        if G.name is None:
            G.name = "candidates"

    if K is not None:
        idx0 = K.coords["sample_0"].values
        idx1 = K.coords["sample_1"].values
        if len(unique(idx0)) < len(idx0) or len(unique(idx1)) < len(idx1):
            msg = "Non-unique sample ids are not allowed in the covariate array"
            msg += " if the sample ids are not equal nor in the same order."
            raise ValueError(msg)
        if K.name is None:
            K.name = "variance-covariance"

    inc_msg = "The provided outcome and {} arrays are sample-wise incompatible."
    try:
        M = M.sel(sample=y.coords["sample"].values)
    except IndexError as e:
        msg = inc_msg.format("covariates")
        raise ValueError(str(e) + "\n\n" + inc_msg)

    if M.name is None:
        M.name = "covariates"

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

    # nsamples = len(y)
    # if nsamples != M.shape[0]:
    #     msg = inc_msg.format("covariates")
    #     raise ValueError(inc_msg)

    # if K is not None:
    #     if nsamples != K.shape[0] or nsamples != K.shape[1]:
    #         inc_msg = inc_msg.format("covariance")
    #         raise ValueError(inc_msg)

    # if G is not None:
    #     if nsamples != G.shape[0]:
    #         inc_msg = inc_msg.format("candidates")
    #         raise ValueError(inc_msg)

    return dict(y=y, M=M, G=G, K=K)


def _dataarray_upcast(x):
    if isinstance(x, (dd.Series, dd.DataFrame)):
        xidx = x.index.compute()
        x = da.asarray(x)
        x = array_shape_reveal(x)
        x0 = xr.DataArray(x)
        x0.coords[x0.dims[0]] = xidx
        if isinstance(x, dd.DataFrame):
            x0.coords[x0.dims[1]] = x.columns
        x = x0

    if not isinstance(x, xr.DataArray):
        x = xr.DataArray(x, encoding={"dtype": "float64"})

    if x.dtype != np.dtype("float64"):
        x = x.astype("float64")

    if x.ndim < 2:
        x = x.expand_dims("dim_1", 1)
    return x


# def _create_index(vals):
#     from xarray import DataArray

#     a = asarray(vals, object)
#     b = asarray(vals, object)
#     return DataArray(a, dims=["sample"], coords={"sample": b})


def _infer_samples_index(arrs):

    samples = arrs[0].coords[arrs[0].dims[0]].values
    for a in arrs[1:]:
        if len(a.coords[a.dims[0]].values) != len(samples):
            break
        if any(a.coords[a.dims[0]].values != samples):
            break
    else:
        return Counter(samples), []

    samples_sets = [Counter(a.coords[a.dims[0]].values) for a in arrs]

    set_intersection = samples_sets[0]
    for ss in samples_sets[1:]:
        set_intersection = set_intersection & ss

    membership_size = [
        asarray([ss[si] for ss in samples_sets], int) for si in set_intersection
    ]

    valid_samples = Counter()
    invalid_samples = []

    for i, k in enumerate(set_intersection.keys()):
        if sum(membership_size[0] > 1) > 1:
            invalid_samples.append(k)
        else:
            valid_samples[k] = set_intersection[k]

    valid_samples = valid_samples
    invalid_samples = invalid_samples
    # if len(valid_samples) > len(invalid_samples):
    #     invalid_samples = invalid_samples.astype(valid_samples.dtype)
    # else:
    #     valid_samples = valid_samples.astype(invalid_samples.dtype)

    return (valid_samples, invalid_samples)


#     ok = _check_samples_arrays_compatibility(arrs)
#     if not ok:
#         msg = "The provided arrays are sample-wise incompatible."
#         msg += " Please, check the number of rows."
#         raise ValueError(msg)

#     if len(arrs) == 0:
#         return asarray([], int)

#     iarrs = [a for a in arrs if _has_sample_index(a)]
#     if len(iarrs) == 0:
#         return asarray(range(arrs[0].shape[0]))

#     index_set = _index_set_intersection(iarrs)

#     return asarray(asarray(_same_order_if_possible(index_set, iarrs)))


# def _check_samples_arrays_compatibility(arrs):
#     """Check index and number of rows.

#     The arrays are compatible if they have the same number of rows or
#     if every array is indexed. If they have the same number of rows but some
#     arrays are indexed and some not, the indexed arrays need to show the same
#     index order as to be compatible to the non-indexed ones.
#     """
#     s = set([len(a) for a in arrs])
#     if len(s) == 0:
#         return True

#     iarrs = [a for a in arrs if _has_sample_index(a)]
#     if len(arrs) == len(iarrs):
#         return True

#     if len(s) == 1:
#         if len(iarrs) == 0:
#             return True
#         # Make sure index have the same order.
#         index = _get_sample_index(iarrs[0])
#         return all([all(index == _get_sample_index(a)) for a in iarrs])

#     return False


def _has_sample_index(x):
    if hasattr(x, "index"):
        return True
    return hasattr(x, "coords") and "sample" == x.dims[0]


def _get_sample_index(x):
    if hasattr(x, "index"):
        return asarray(x.index.values)
    return asarray(x.coords["sample"].values)


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
    return asarray([], int)
