from collections import Counter

import numpy as np
from numpy import array_equal, asarray, unique

from ._dask import array_shape_reveal


def normalise_dataset(y, M=None, G=None, K=None):
    r"""Convert data types to DataArray.

    This is a fundamental function for :mod:`limix` as it standardise outcome,
    covariates, candidates, and kinship arrays into :class:`xarray.DataArray` data type.
    Data arrays are :mod:`numpy`/:mod:`dask` arrays with indexed coordinates,
    therefore generalising data frames from :mod:`pandas`. It allows for lazy loading of
    data via dask arrays. It also supports arrays with different dimensionality and
    types, mixture of indexed and non-indexed arrays, and repeated sample labels.

    Examples
    --------

    .. doctest::

        >>> import pytest
        >>> from numpy.random import RandomState
        >>> from pandas import DataFrame
        >>> from xarray import DataArray
        >>> from limix._dataset import normalise_dataset
        >>>
        >>> random = RandomState(0)
        >>>
        >>> y = random.randn(4)
        >>> y = DataFrame(y, index=["sample0", "sample0", "sample1", "sample2"])
        >>>
        >>> G = random.randn(5, 6)
        >>>
        >>> data = normalise_dataset(y, G=G)
        >>> print(data["y"])
        <xarray.DataArray 'outcome' (sample: 4, trait: 1)>
        array([[1.764052],
               [0.400157],
               [0.978738],
               [2.240893]])
        Coordinates:
          * sample   (sample) object 'sample0' 'sample0' 'sample1' 'sample2'
          * trait    (trait) int64 0
        >>> print(data["G"])
        <xarray.DataArray 'candidates' (sample: 4, candidate: 6)>
        array([[ 1.867558, -0.977278,  0.950088, -0.151357, -0.103219,  0.410599],
               [ 0.144044,  1.454274,  0.761038,  0.121675,  0.443863,  0.333674],
               [ 1.494079, -0.205158,  0.313068, -0.854096, -2.55299 ,  0.653619],
               [ 0.864436, -0.742165,  2.269755, -1.454366,  0.045759, -0.187184]])
        Coordinates:
          * sample   (sample) object 'sample0' 'sample0' 'sample1' 'sample2'
        Dimensions without coordinates: candidate
        >>> K = random.randn(3, 3)
        >>> K = K.dot(K.T)
        >>> K = DataArray(K)
        >>> K.coords["dim_0"] = ["sample0", "sample1", "sample2"]
        >>> K.coords["dim_1"] = ["sample0", "sample1", "sample2"]
        >>>
        >>> data = normalise_dataset(y, K=K)
        >>> print(data["y"])
        <xarray.DataArray 'outcome' (sample: 4, trait: 1)>
        array([[1.764052],
               [0.400157],
               [0.978738],
               [2.240893]])
        Coordinates:
          * sample   (sample) object 'sample0' 'sample0' 'sample1' 'sample2'
          * trait    (trait) int64 0
        >>> print(data["K"])
        <xarray.DataArray 'variance-covariance' (sample_0: 4, sample_1: 4)>
        array([[ 1.659103,  1.659103, -0.850801, -1.956422],
               [ 1.659103,  1.659103, -0.850801, -1.956422],
               [-0.850801, -0.850801,  1.687126, -0.194938],
               [-1.956422, -1.956422, -0.194938,  6.027272]])
        Coordinates:
          * sample_0  (sample_0) <U7 'sample0' 'sample0' 'sample1' 'sample2'
          * sample_1  (sample_1) <U7 'sample0' 'sample0' 'sample1' 'sample2'
        >>> with pytest.raises(ValueError):
        ...     normalise_dataset(y, G=G, K=K)
    """

    y = _rename_dims(_dataarray_upcast(y), "sample", "trait")
    M = _rename_dims(_dataarray_upcast(M), "sample", "covariate")
    G = _rename_dims(_dataarray_upcast(G), "sample", "candidate")
    K = _rename_dims(_dataarray_upcast(K), "sample_0", "sample_1")

    data = {"y": y, "M": M, "G": G, "K": K}
    dim_name = {
        "y": ["sample"],
        "M": ["sample"],
        "G": ["sample"],
        "K": ["sample_0", "sample_1"],
    }
    arrname = {
        "y": "outcome",
        "M": "covariates",
        "G": "candidates",
        "K": "variance-covariance",
    }

    data.update(_assign_index_to_nonindexed(_fout(data), dim_name))

    valid_samples = _infer_samples_index(_fout(data), dim_name)
    for k, v in data.items():
        for dn in dim_name[k]:
            data[k] = _assign_coords(v, dn, valid_samples)

    for k, v in _fout(data).items():
        if v.name is None:
            v.name = arrname[k]

    n = len(data["y"].coords["sample"])
    o = [v.coords[dn].size == n for k, v in _fout(data).items() for dn in dim_name[k]]
    if all(o):
        return data

    _check_uniqueness(data, dim_name, arrname)
    _check_sample_compatibility(data, dim_name, arrname)

    return data


def _assign_index_to_nonindexed(data, dim_name):
    indexed = {}
    nonindexed = {}

    for k, v in data.items():
        for dn in dim_name[k]:
            if dn not in v.coords:
                nonindexed[k] = v
            else:
                indexed[k] = v

    if len(nonindexed) == 0:
        return data

    if len(indexed) == 0:
        for k, v in data.items():
            for dn in dim_name[k]:
                n = len(v.coords[dn])
                v.coords[dn] = ["sample{}".format(j) for j in range(n)]
        return data

    k = next(iter(indexed.keys()))
    index = indexed[k].coords[dim_name[k][0]]
    for k, v in indexed.items():
        for dn in dim_name[k]:
            if not array_equal(index, v.coords[dn]):
                msg = "Please, check the provided sample labels in your arrays."
                msg = " There are some inconsistences in them."
                raise ValueError(msg)

    n = min(v.coords[dn].size for (k, v) in data.items() for dn in dim_name[k])
    index = index[:n]

    for k, v in data.items():
        for dn in dim_name[k]:
            data[k] = _take(v, dn, n)

    nonindexed = {}
    for k, v in data.items():
        for dn in dim_name[k]:
            if dn not in v.coords:
                nonindexed[k] = v

    for k, v in nonindexed.items():
        for dn in dim_name[k]:
            v.coords[dn] = index.values

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


def _infer_samples_index(data, dim_name):

    k, v = next(iter(data.items()))
    samples = v.coords[dim_name[k][0]].values
    for k, v in data.items():
        for dn in dim_name[k]:
            if not array_equal(v.coords[dn].values, samples):
                break
        else:
            continue
        break
    else:
        return Counter(samples)

    samples_sets = [
        Counter(v.coords[dn].values) for k, v in data.items() for dn in dim_name[k]
    ]

    set_intersection = samples_sets[0]
    for ss in samples_sets[1:]:
        set_intersection = set_intersection & ss

    membership_size = [
        asarray([ss[si] for ss in samples_sets], int) for si in set_intersection
    ]

    valid_samples = Counter()

    for i, k in enumerate(set_intersection.keys()):
        if sum(membership_size[0] > 1) <= 1:
            valid_samples[k] = set_intersection[k]

    return valid_samples


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


def _fout(data):
    return {k: v for k, v in data.items() if v is not None}


def _take(x, dim_name, n):
    i = 0
    sl = [slice(None)] * x.ndim
    while i < x.ndim and x.dims[i] != dim_name:
        i += 1
    sl[i] = slice(0, n)
    return x[tuple(sl)]


def _check_uniqueness(data, dim_name, arrname):
    msg = "Non-unique sample ids are not allowed in the {} array"
    msg += " if the sample ids are not equal nor in the same order."

    for k, v in _fout(data).items():
        if k == "y":
            continue
        for dn in dim_name[k]:
            idx = v.coords[dn].values
            if len(unique(idx)) < len(idx):
                raise ValueError(msg.format(arrname[k]))


def _check_sample_compatibility(data, dim_name, arrname):
    inc_msg = "The provided outcome and {} arrays are sample-wise incompatible."

    for k in _fout(data).keys():
        if k == "y":
            continue
        for dn in dim_name[k]:
            try:
                data[k] = data[k].sel(**{dn: data["y"].coords["sample"].values})
            except IndexError as e:
                raise ValueError(str(e) + "\n\n" + inc_msg.format(arrname[k]))
