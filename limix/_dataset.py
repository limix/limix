from collections import Counter

import numpy as np
from numpy import array_equal, asarray, unique

from ._dask import array_shape_reveal


def normalise_dataset(y, M=None, G=None, K=None):

    y = _dataarray_upcast(y)
    y = y.rename({y.dims[0]: "sample"})
    y = y.rename({y.dims[1]: "trait"})

    if M is not None:
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

    data = {"y": y, "M": M, "G": G, "K": K}
    datas = {k: (data[k], 0) for k in ["y", "M", "G"] if data[k] is not None}
    if data["K"] is not None:
        datas["K0"] = (K, 0)
        datas["K1"] = (K, 1)

    datas = _assign_samples_index_to_nonindexed(datas)
    valid_samples, invalid_samples = _infer_samples_index(list(datas.values()))

    y = datas["y"][0]

    if "M" in datas:
        M = datas["M"][0]

    if "G" in datas:
        G = datas["G"][0]

    if "K0" in datas:
        K = datas["K0"][0]

    if "sample" not in y.coords:
        y = y.assign_coords(sample=list(valid_samples))
    else:
        y = y.loc[{"sample": y.get_index("sample").isin(valid_samples)}]

    if M is not None:
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

    if M is not None:
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
    if M is not None:
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

    return dict(y=y, M=M, G=G, K=K)


def _assign_samples_index_to_nonindexed(datas):
    indexed = []
    nonindexed = []

    for a in datas.values():
        v, i = a
        if v.dims[i] not in v.coords:
            nonindexed.append(a)
        else:
            indexed.append(a)

    if len(nonindexed) == 0:
        return datas

    if len(indexed) == 0:
        for a in datas.values():
            v, i = a
            n = v.shape[i]
            v.coords[v.dims[i]] = ["sample{}".format(j) for j in range(n)]
        return datas

    v, i = indexed[0]
    index = v.coords[v.dims[i]]
    for a in indexed[1:]:
        v, i = a
        if not array_equal(index, v.coords[v.dims[i]]):
            msg = "Please, check the provided sample labels in your arrays."
            msg = " There are some inconsistences in them."
            raise ValueError(msg)

    # n = min(v.coords[v.dims[i]].size for (v, i) in arrs)
    n = min(v.coords[v.dims[i]].size for (v, i) in datas.values())
    index = index[:n]

    for (k, (v, i)) in datas.items():
        sl = [slice(None)] * v.ndim
        sl[i] = slice(0, n)
        v = v[tuple(sl)]
        datas[k] = (v, i)

    indexed = []
    nonindexed = []
    for (v, i) in datas.values():
        if v.dims[i] not in v.coords:
            nonindexed.append((v, i))
        else:
            indexed.append((v, i))

    for (v, i) in nonindexed:
        v.coords[v.dims[i]] = index.values

    return datas


def _dataarray_upcast(x):
    import dask.dataframe as dd
    import dask.array as da
    import xarray as xr

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


def _infer_samples_index(arrs):

    v, i = arrs[0]
    samples = v.coords[v.dims[i]].values
    for a in arrs[1:]:
        v, i = a
        if not array_equal(v.coords[v.dims[i]].values, samples):
            break
    else:
        return Counter(samples), []

    samples_sets = [Counter(a[0].coords[a[0].dims[a[1]]].values) for a in arrs]

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

    return (valid_samples, invalid_samples)


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
