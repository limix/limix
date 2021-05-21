from collections import Counter

from .._bits.deco import return_none_if_none
from .._bits.xarray import set_coord
from ._asarray import asarray as _asarray
from ._conf import CONF

set_coord = return_none_if_none(set_coord)
_asarray = return_none_if_none(_asarray)


def conform_dataset(y, M=None, G=None, K=None):
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
        array([[1.76405235],
               [0.40015721],
               [0.97873798],
               [2.24089320]])
        Coordinates:
          * sample   (sample) object 'sample0' 'sample0' 'sample1' 'sample2'
          * trait    (trait) int64 0
        >>> print(data["G"])
        <xarray.DataArray 'genotype' (sample: 4, candidate: 6)>
        array([[ 1.86755799, -0.97727788,  0.95008842, -0.15135721, -0.10321885,
                0.41059850],
              [ 0.14404357,  1.45427351,  0.76103773,  0.12167502,  0.44386323,
                0.33367433],
              [ 1.49407907, -0.20515826,  0.31306770, -0.85409574, -2.55298982,
                0.65361860],
              [ 0.86443620, -0.74216502,  2.26975462, -1.45436567,  0.04575852,
                -0.18718385]])
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
        array([[1.76405235],
            [0.40015721],
            [0.97873798],
            [2.24089320]])
        Coordinates:
        * sample   (sample) object 'sample0' 'sample0' 'sample1' 'sample2'
        * trait    (trait) int64 0
        >>> print(data["K"])
        <xarray.DataArray 'covariance' (sample_0: 4, sample_1: 4)>
        array([[ 1.65910302,  1.65910302, -0.85080096, -1.95642234],
               [ 1.65910302,  1.65910302, -0.85080096, -1.95642234],
               [-0.85080096, -0.85080096,  1.68712632, -0.19493776],
               [-1.95642234, -1.95642234, -0.19493776,  6.02727223]])
        Coordinates:
        * sample_0  (sample_0) object 'sample0' 'sample0' 'sample1' 'sample2'
        * sample_1  (sample_1) object 'sample0' 'sample0' 'sample1' 'sample2'
        >>> with pytest.raises(ValueError):
        ...     conform_dataset(y, G=G, K=K)
    """
    y = _asarray(y, "trait", CONF["data_dims"]["trait"])
    M = _asarray(M, "covariate", CONF["data_dims"]["covariate"])
    G = _asarray(G, "genotype", CONF["data_dims"]["genotype"])
    K = _asarray(K, "covariance", CONF["data_dims"]["covariance"])

    data = {"y": y, "M": M, "G": G, "K": K}
    data = {k: v for k, v in data.items() if v is not None}

    sample_dims = [
        t
        for t in [
            ("y", "sample"),
            ("M", "sample"),
            ("G", "sample"),
            ("K", "sample_0"),
            ("K", "sample_1"),
        ]
        if t[0] in data
    ]

    data = _fix_samples(data, sample_dims)
    for n in data.keys():
        data[n].name = CONF["varname_to_target"][n]

    nsamples = len(data["y"].coords["sample"])
    same_size = all(data[n].coords[d].size == nsamples for n, d in sample_dims)

    data = _fix_covariates(data, same_size)

    # We accept non-unique samples when all sample sizes are equal.
    # In the other cases, we check for uniqueness.
    if not same_size:
        _check_uniqueness(data, sample_dims)
    _match_samples(data, sample_dims)

    return {k: data.get(k, None) for k in ["y", "M", "G", "K"]}


def _default_covariates(samples):
    from numpy import ones, asarray
    from xarray import DataArray

    M = ones((samples.size, 1))
    M = DataArray(
        M,
        # encoding={"dtype": "float64"},
        dims=["sample", "covariate"],
        coords={"sample": samples, "covariate": asarray(["offset"], dtype=object)},
    )

    return M


def _default_sample_coords(n):
    return ["sample{}".format(j) for j in range(n)]


def _fix_samples(data, sample_dims):
    from .._bits.xarray import take
    from numpy import array_equal

    samples_list = [
        data[n].coords[d].values for n, d in sample_dims if d in data[n].coords
    ]
    nmin_samples = min(data[n].coords[d].size for n, d in sample_dims)

    if not samples_list:
        data["y"] = take(data["y"], slice(0, nmin_samples), "sample")
        data["y"].coords["sample"] = _default_sample_coords(data["y"].sample.size)
        samples_list.append(data["y"].coords["sample"].values)

    samples = samples_list[0]

    if len(samples_list) < len(sample_dims):
        if not all([array_equal(s, samples) for s in samples_list]):
            raise ValueError(
                "Please, check the provided sample labels in your arrays."
                + " There are some inconsistences between them."
            )

        for n, d in sample_dims:
            data[n] = take(data[n], slice(0, nmin_samples), d)
            data[n].coords[d] = samples[:nmin_samples]

    samples_list = [
        data[n].coords[d].values for n, d in sample_dims if d in data[n].coords
    ]

    valid_samples = _infer_samples_index(samples_list)
    for n, d in sample_dims:
        data[n] = set_coord(data[n], d, valid_samples)

    return data


def _infer_samples_index(samples_list):
    r""" Infer a list of sample labels that is compatible to the provided data.

    It uses :class:`collections.Counter` to count the number of repeated sample
    labels, and to provide set (bag) intersection that handles repeated elements.
    """
    from numpy import array_equal, asarray

    samples = samples_list[0]
    if all(array_equal(s, samples) for s in samples_list):
        return Counter(samples)

    samples_sets = [Counter(s) for s in samples_list]

    set_intersection = samples_sets[0]
    for ss in samples_sets[1:]:
        set_intersection = set_intersection & ss

    membership_size = [
        asarray([ss[si] for ss in samples_sets], int) for si in set_intersection
    ]

    valid_samples = Counter()

    for k in set_intersection.keys():
        if sum(membership_size[0] > 1) <= 1:
            valid_samples[k] = set_intersection[k]

    return valid_samples


def _fix_covariates(data, samples_same_size):
    from pandas import unique

    y = data["y"]

    # We accept non-unique samples when all sample sizes are equal.
    if samples_same_size:
        if "M" not in data:
            data["M"] = _default_covariates(y.sample.values)
    else:
        if "M" not in data:
            data["M"] = _default_covariates(unique(y.sample.values))

    return data


def _check_uniqueness(data, dims):
    from numpy import unique

    msg = "Non-unique sample ids are not allowed in the {} array"
    msg += " if the sample ids are not equal nor in the same order."

    for n, d in dims:
        if n == "y":
            continue
        idx = data[n].coords[d].values
        if len(unique(idx)) < len(idx):
            raise ValueError(msg.format(data[n].name))


def _match_samples(data, dims):
    from pandas.core.indexes.api import InvalidIndexError

    inc_msg = "The provided trait and {} arrays are sample-wise incompatible."

    for n, d in dims:
        if n == "y":
            continue
        try:
            try:
                data[n] = data[n].sel(**{d: data["y"].coords["sample"].values})
            except InvalidIndexError:
                pass
        except IndexError as e:
            raise ValueError(str(e) + "\n\n" + inc_msg.format(data[n].name))

    return data
