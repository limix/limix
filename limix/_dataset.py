from collections import Counter

import numpy as np
from numpy import array_equal, asarray, unique

from ._dask import array_shape_reveal


def normalise_dataset(y, M=None, G=None, K=None):

    y = _rename_dims(_dataarray_upcast(y), "sample", "trait")
    M = _rename_dims(_dataarray_upcast(M), "sample", "covariate")
    G = _rename_dims(_dataarray_upcast(G), "sample", "candidate")
    K = _rename_dims(_dataarray_upcast(K), "sample_0", "sample_1")

    data = {"y": y, "M": M, "G": G, "K0": K, "K1": K}
    dim = {"y": 0, "M": 0, "G": 0, "K0": 0, "K1": 1}

    data = _assign_index_to_nonindexed(
        {k: v for (k, v) in data.items() if v is not None}, dim
    )
    for k in dim.keys():
        if k not in data:
            data[k] = None

    valid_samples, invalid_samples = _infer_samples_index(
        {k: v for (k, v) in data.items() if v is not None}, dim
    )

    data["y"] = _assign_coords(data["y"], "sample", valid_samples)
    data["M"] = _assign_coords(data["M"], "sample", valid_samples)
    data["G"] = _assign_coords(data["G"], "sample", valid_samples)
    data["K0"] = _assign_coords(data["K0"], "sample_0", valid_samples)
    data["K0"] = _assign_coords(data["K0"], "sample_1", valid_samples)
    data["K1"] = data["K0"]

    n = len(data["y"].coords["sample"])
    if data["y"].name is None:
        data["y"].name = "outcome"

    names = {
        "M": "covariates",
        "G": "candidates",
        "K0": "variance-covariance",
        "K1": "variance-covariance",
    }
    for k, name in names.items():
        v = data[k]
        if v is None:
            continue
        if v.name is None:
            v.name = name

    o = [len(v.coords[v.dims[dim[k]]]) == n for (k, v) in data.items() if v is not None]

    if all(o):
        del data["K1"]
        data["K"] = data["K0"]
        del data["K0"]
        return data

    msg = "Non-unique sample ids are not allowed in the {} array"
    msg += " if the sample ids are not equal nor in the same order."

    for k, name in names.items():
        v = data[k]
        if v is None:
            continue

        idx = v.coords[v.dims[dim[k]]].values

        if len(unique(idx)) < len(idx):
            raise ValueError(msg.format(name))

    inc_msg = "The provided outcome and {} arrays are sample-wise incompatible."

    dimn = {"sample": ["y", "M", "G"], "sample_0": ["K"], "sample_1": ["K"]}
    del data["K1"]
    data["K"] = data["K0"]
    del data["K0"]

    for lab, li in dimn.items():
        for k in li:
            if k == "y":
                continue
            v = data[k]
            if v is None:
                continue
            try:
                data[k] = v.sel(**{lab: data["y"].coords["sample"].values})
            except IndexError as e:
                raise ValueError(str(e) + "\n\n" + inc_msg.format(name))

    return data


def _assign_index_to_nonindexed(data, dim):
    indexed = {}
    nonindexed = {}

    for k, v in data.items():
        if v.dims[dim[k]] not in v.coords:
            nonindexed[k] = v
        else:
            indexed[k] = v

    if len(nonindexed) == 0:
        return data

    if len(indexed) == 0:
        for k, v in data.items():
            n = v.shape[dim[k]]
            v.coords[v.dims[dim[k]]] = ["sample{}".format(j) for j in range(n)]
        return data

    k, v = next(iter(indexed.items()))
    index = v.coords[v.dims[dim[k]]]
    for k, v in indexed.items():
        if not array_equal(index, v.coords[v.dims[dim[k]]]):
            msg = "Please, check the provided sample labels in your arrays."
            msg = " There are some inconsistences in them."
            raise ValueError(msg)

    n = min(v.coords[v.dims[dim[k]]].size for (k, v) in data.items())
    index = index[:n]

    for k, v in data.items():
        sl = [slice(None)] * v.ndim
        sl[dim[k]] = slice(0, n)
        v = v[tuple(sl)]
        data[k] = v

    nonindexed = {}
    for k, v in data.items():
        if v.dims[dim[k]] not in v.coords:
            nonindexed[k] = v

    for k, v in nonindexed.items():
        v.coords[v.dims[dim[k]]] = index.values

    return data


def _rename_dims(x, dim_0, dim_1):
    if x is None:
        return None
    x = x.rename({x.dims[0]: dim_0})
    x = x.rename({x.dims[1]: dim_1})
    return x


def _dataarray_upcast(x):
    import dask.dataframe as dd
    import dask.array as da
    import xarray as xr

    if x is None:
        return None

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


def _assign_coords(x, dim_name, samples):
    if x is None:
        return None
    if dim_name not in x.coords:
        x = x.assign_coords(**{dim_name: list(samples)})
    else:
        x = x.loc[{dim_name: x.get_index(dim_name).isin(samples)}]
    return x


def _infer_samples_index(data, dim):

    k, v = next(iter(data.items()))
    samples = v.coords[v.dims[dim[k]]].values
    for k, v in data.items():
        if not array_equal(v.coords[v.dims[dim[k]]].values, samples):
            break
    else:
        return Counter(samples), []

    samples_sets = [Counter(v.coords[v.dims[dim[k]]].values) for (k, v) in data.items()]

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
