from __future__ import unicode_literals
from collections import Counter
from functools import wraps

from numpy import array_equal, asarray, unique, dtype

from .._bits.dask import array_shape_reveal
from .._bits.xarray import set_coord
from ._dataarray import fix_dim_order_if_hinted, rename_dims
from .conf import is_data_name, is_short_data_name, data_name


def conform_dataset(y, M=None, G=None, K=None, X=None):
    r""" Convert data types to DataArray.

    This is a fundamental function for :mod:`limix` as it standardise outcome,
    covariates, genotype, and kinship arrays into :class:`xarray.DataArray` data type.
    Data arrays are :mod:`numpy`/:mod:`dask` arrays with indexed coordinates,
    therefore generalising data frames from :mod:`pandas`. It allows for lazy loading of
    data via dask arrays. It also supports arrays with different dimensionality and
    types, mixture of indexed and non-indexed arrays, and repeated sample labels.

    Examples
    --------

    .. doctest::

        >>> from __future__ import unicode_literals
        >>> import pytest
        >>> from numpy.random import RandomState
        >>> from pandas import DataFrame
        >>> from xarray import DataArray
        >>> from limix._data import conform_dataset
        >>>
        >>> random = RandomState(0)
        >>>
        >>> y = random.randn(4)
        >>> y = DataFrame(y, index=["sample0", "sample0", "sample1", "sample2"])
        >>>
        >>> G = random.randn(5, 6)
        >>>
        >>> data = conform_dataset(y, G=G)
        >>> print(data["y"])
        <xarray.DataArray 'trait' (sample: 4, trait: 1)>
        array([[1.764052],
               [0.400157],
               [0.978738],
               [2.240893]])
        Coordinates:
          * sample   (sample) object 'sample0' 'sample0' 'sample1' 'sample2'
          * trait    (trait) int64 0
        >>> print(data["G"])
        <xarray.DataArray 'genotype' (sample: 4, candidate: 6)>
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
        >>> data = conform_dataset(y, K=K)
        >>> print(data["y"])
        <xarray.DataArray 'trait' (sample: 4, trait: 1)>
        array([[1.764052],
               [0.400157],
               [0.978738],
               [2.240893]])
        Coordinates:
          * sample   (sample) object 'sample0' 'sample0' 'sample1' 'sample2'
          * trait    (trait) int64 0
        >>> print(data["K"])
        <xarray.DataArray 'covariance' (sample_0: 4, sample_1: 4)>
        array([[ 1.659103,  1.659103, -0.850801, -1.956422],
               [ 1.659103,  1.659103, -0.850801, -1.956422],
               [-0.850801, -0.850801,  1.687126, -0.194938],
               [-1.956422, -1.956422, -0.194938,  6.027272]])
        Coordinates:
          * sample_0  (sample_0) object 'sample0' 'sample0' 'sample1' 'sample2'
          * sample_1  (sample_1) object 'sample0' 'sample0' 'sample1' 'sample2'
        >>> with pytest.raises(ValueError):
        ...     conform_dataset(y, G=G, K=K)
    """
    from pandas import unique

    if X is None:
        X = []

    y = rename_dims(fix_dim_order_if_hinted(to_dataarray(y)), ["sample", "trait"])
    M = rename_dims(fix_dim_order_if_hinted(to_dataarray(M)), ["sample", "covariate"])
    G = rename_dims(fix_dim_order_if_hinted(to_dataarray(G)), ["sample", "candidate"])
    K = rename_dims(fix_dim_order_if_hinted(to_dataarray(K)), ["sample_0", "sample_1"])

    X = [
        (rename_dims(fix_dim_order_if_hinted(to_dataarray(x)), ["sample", n1]), n0, n1)
        for x, n0, n1 in X
    ]

    data = {"y": y, "M": M, "G": G, "K": K}
    data.update({"X{}".format(i): x[0] for i, x in enumerate(X)})

    dim_name = {
        "y": ["sample"],
        "M": ["sample"],
        "G": ["sample"],
        "K": ["sample_0", "sample_1"],
    }
    dim_name.update({"X{}".format(i): ["sample"] for i in range(len(X))})

    samples_list = []
    missing_samples = False
    for k, v in _fout(data).items():
        for dn in dim_name[k]:
            if dn in v.coords:
                samples_list.append(v.coords[dn])
            else:
                missing_samples = True

    if len(samples_list) == 0:
        data.update(_assign_coords(_fout(data), dim_name, create_default_sample_coords))
    elif missing_samples:
        # samples = min(samples_list, key=lambda s: s.size)
        samples = samples_list[0].values
        ok = all([array_equal(s, samples) for s in samples_list])
        if not ok:
            msg = "Please, check the provided sample labels in your arrays."
            msg += " There are some inconsistences between them."
            raise ValueError(msg)
        n = min(v.coords[dn].size for (k, v) in _fout(data).items() for dn in dim_name[k])
        samples = samples[:n]
        data.update(_slice_coordinates(_fout(data), dim_name, len(samples)))
        # data.update(_set_coordinates(_fout(data), dim_name, len(samples)))
        data.update(_assign_coords(_fout(data), dim_name, lambda _: samples))
        # data.update(_assign_index_to_nonindexed(_fout(data), dim_name))

    valid_samples = _infer_samples_index(_fout(data), dim_name)

    for k, v in data.items():
        for dn in dim_name[k]:
            data[k] = set_coord(v, dn, valid_samples)

    for k, v in _fout(data).items():
        if v.name is None:
            if is_short_data_name(k) or is_data_name(k):
                v.name = data_name(k)
            else:
                v.name = k

    n = len(data["y"].coords["sample"])
    o = [v.coords[dn].size == n for k, v in _fout(data).items() for dn in dim_name[k]]

    if all(o):
        if data["M"] is None:
            data["M"] = create_default_covariates(y.sample.values)
        return data

    if data["M"] is None:
        data["M"] = create_default_covariates(unique(y.sample.values))

    _check_uniqueness(data, dim_name)
    _check_sample_compatibility(data, dim_name)

    return data


def to_dataarray(x):
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

    if x.dtype != dtype("float64"):
        x = x.astype("float64")

    if x.ndim < 2:
        x = x.expand_dims("dim_1", 1)

    for dim in x.dims:
        if x.coords[dim].dtype.kind in {"U", "S"}:
            x.coords[dim] = x.coords[dim].values.astype(object)
    return x


def create_default_covariates(samples):
    from numpy import ones, asarray
    from xarray import DataArray

    M = ones((samples.size, 1))
    M = DataArray(
        M,
        encoding={"dtype": "float64"},
        dims=["sample", "covariate"],
        coords={"sample": samples, "covariate": asarray(["offset"], dtype=object)},
    )

    return M


def create_default_sample_coords(n):
    return ["sample{}".format(j) for j in range(n)]


def _assign_index_to_nonindexed(data, dims):
    from .._bits.xarray import take

    present, absent = _split_coordinates(data, dims)
    if len(absent) == 0:
        return data

    if not _check_equal_coords(present, dims):
        msg = "Please, check the provided sample labels in your arrays."
        msg += " There are some inconsistences between them."
        raise ValueError(msg)

    # The reference index is the smallest one we can find.
    n = min(v.coords[dn].size for (k, v) in data.items() for dn in dims[k])
    k = next(iter(present.keys()))
    index = present[k].coords[dims[k][0]][:n]

    # Take only the first n indices.
    # The rest will be ignored.
    for k, v in data.items():
        for dn in dims[k]:
            data[k] = take(v, slice(0, n), dn)

    _, absent = _split_coordinates(data, dims)

    for k, v in absent.items():
        for dn in dims[k]:
            v.coords[dn] = index.values

    return data


def _get_shortest_coords(data, dims):
    n = min(v.coords[dn].size for (k, v) in data.items() for dn in dims[k])
    k = next(iter(data.keys()))
    return data[k].coords[dims[k][0]][:n]


def _split_coordinates(data, dims):
    r""" Split the mentioned coordinates into absent or present. """
    present = {}
    absent = {}
    for name, c in data.items():
        for dim in dims[name]:
            if dim not in c.coords:
                absent[name] = c
            else:
                present[name] = c
    return present, absent


def _slice_coordinates(data, dims, n):
    from .._bits.xarray import take

    for name, c in data.items():
        for dim in dims[name]:
            data[name] = take(c, slice(0, n), dim)
    return data


def _any_coords(data, dims):
    r""" Check if `data` has any data array with the specified dimension name. """
    present, _ = _split_coordinates(data, dims)
    return len(present) > 0


def _all_coords(data, dims):
    r""" Check if every data array has its specified dimension name. """
    _, absent = _split_coordinates(data, dims)
    return len(absent) == 0


def _assign_coords(data, dims, create_coords):
    r""" Assign coordinate values for the specified dimension name. """
    for name, c in data.items():
        for dim in dims[name]:
            n = len(c.coords[dim])
            c.coords[dim] = create_coords(n)
    return data


def _check_equal_coords(data, dims):
    name = next(iter(data.keys()))
    coords = data[name].coords[dims[name][0]]

    for name, c in data.items():
        for dim in dims[name]:
            if not array_equal(coords, c.coords[dim]):
                return False
    return True


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


def _fout(data):
    return {k: v for k, v in data.items() if v is not None}


def _check_uniqueness(data, dim_name):
    msg = "Non-unique sample ids are not allowed in the {} array"
    msg += " if the sample ids are not equal nor in the same order."

    for k, v in _fout(data).items():
        if k == "y":
            continue
        for dn in dim_name[k]:
            idx = v.coords[dn].values
            if len(unique(idx)) < len(idx):
                raise ValueError(msg.format(v.name))


def _check_sample_compatibility(data, dim_name):
    inc_msg = "The provided trait and {} arrays are sample-wise incompatible."

    for k, v in _fout(data).items():
        if k == "y":
            continue
        for dn in dim_name[k]:
            try:
                data[k] = data[k].sel(**{dn: data["y"].coords["sample"].values})
            except IndexError as e:
                raise ValueError(str(e) + "\n\n" + inc_msg.format(v.name))


def _return_none(f):
    @wraps(f)
    def new_func(x, *args, **kwargs):
        if x is None:
            return None
        return f(x, *args, **kwargs)

    return new_func


rename_dims = _return_none(rename_dims)
to_dataarray = _return_none(to_dataarray)
fix_dim_order_if_hinted = _return_none(fix_dim_order_if_hinted)
set_coord = _return_none(set_coord)
