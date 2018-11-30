from numpy import array, asarray, dtype
from numpy.random import RandomState
from numpy.testing import assert_, assert_array_equal, assert_equal

from limix._data import conform_dataset
from pandas import DataFrame, Series
from xarray import DataArray


def test_dataset_conform_dataset():
    y = array([-1.2, 3.4, 0.1])
    samples = ["sample{}".format(i) for i in range(len(y))]

    y = DataFrame(data=y, index=samples)

    random = RandomState(0)

    K = random.randn(3, 4)
    K = K.dot(K.T)
    K = DataFrame(data=K, index=samples, columns=samples)

    M = random.randn(3, 2)
    M = DataFrame(data=M, index=samples)

    G = random.randn(2, 4)
    G = DataFrame(data=G, index=samples[:2])

    data = conform_dataset(y, M=M, K=K)

    assert_array_equal(y.values, data["y"].values)

    y = array([-1.2, 3.4, 0.1, 0.1, 0.0, -0.2])

    data = conform_dataset(DataFrame(data=y, index=samples + samples), M=M, G=G, K=K)

    assert_equal(data["y"].shape, (4, 1))
    assert_equal(data["M"].shape, (4, 2))
    assert_equal(data["G"].shape, (4, 4))
    assert_equal(data["K"].shape, (4, 4))

    samples = ["sample0", "sample1", "sample0", "sample1"]
    assert_array_equal(data["y"].sample, samples)
    assert_array_equal(data["M"].sample, samples)
    assert_array_equal(data["G"].sample, samples)
    assert_array_equal(data["K"].sample_0, samples)
    assert_array_equal(data["K"].sample_1, samples)

    assert_array_equal(data["M"].covariate, [0, 1])
    assert_array_equal(data["G"].candidate, [0, 1, 2, 3])


def test_dataset_pandas_xarray_dask():
    import numpy as np
    import dask.array as da
    import dask.dataframe as dd
    import pandas as pd
    from limix._data._conform import to_dataarray

    x = []

    x.append([1.0, 2.0, 3.0])
    x.append(np.asarray([1.0, 2.0, 3.0]))
    x.append(np.asarray([[1.0], [2.0], [3.0]]))
    x.append(np.asarray([[1], [2], [3]], dtype=int))
    x.append(da.from_array(x[0], 2))
    x.append(da.from_array(x[1], 2))
    x.append(da.from_array(x[2], 2))
    x.append(da.from_array(x[3], 2))

    n = len(x)
    for i in range(n):
        if isinstance(x[i], da.Array):
            tmp = np.asarray(x[i])
            if tmp.ndim == 2:
                tmp = tmp.ravel()
                x.append(dd.from_array(tmp))
            else:
                x.append(dd.from_array(x[i]))
        else:
            tmp = np.asarray(x[i])
            if tmp.ndim == 2:
                tmp = tmp.ravel()
                x.append(pd.Series(tmp))
            else:
                x.append(pd.Series(x[i]))

    for i in range(n):
        if isinstance(x[i], da.Array):
            x.append(dd.from_array(x[i]))
        elif isinstance(x[i], np.ndarray):
            x.append(pd.DataFrame(x[i]))

    n = len(x)

    for i in range(n):
        x.append(DataArray(x[i]))
        x.append(x[-1].chunk(2))

    print()
    for xi in x:
        y = to_dataarray(xi)
        assert_equal(y.dtype, dtype("float64"))
        assert_array_equal(y.shape, (3, 1))
        assert_(isinstance(y, DataArray))
        if isinstance(xi, Series):
            assert_array_equal(list(xi.index), list(y.coords["dim_0"].values))
        if isinstance(xi, DataFrame):
            assert_array_equal(list(xi.columns), list(y.coords["dim_1"].values))

        is_dask = (
            hasattr(xi, "chunks")
            and xi.chunks is not None
            or hasattr(xi, "values")
            and hasattr(xi, "values")
            and hasattr(xi.values, "chunks")
            and xi.values.chunks is not None
        )

        assert_equal(is_dask, y.chunks is not None)
        assert_array_equal(asarray(xi).ravel(), asarray(y).ravel())


def test_dataset_different_size():
    random = RandomState(0)
    n0 = 5
    n1 = 3
    y = random.randn(n0)
    samples = ["sample{}".format(i) for i in range(len(y))]
    y = DataFrame(data=y, index=samples)

    G = random.randn(n1, 10)

    data = conform_dataset(y, G=G)

    assert_array_equal(data["y"].values, y[:n1])
    assert_array_equal(data["G"].values, G[:n1, :])

    n0 = 3
    n1 = 5
    y = random.randn(n0)
    samples = ["sample{}".format(i) for i in range(len(y))]
    y = DataFrame(data=y, index=samples)

    G = random.randn(n1, 10)

    data = conform_dataset(y, G=G)

    assert_array_equal(data["y"].values, y[:n0])
    assert_array_equal(data["G"].values, G[:n0, :])


def test_dataset_underline_prefix():

    data = {
        "coords": {
            "trait": {"data": "gene1", "dims": (), "attrs": {}},
            "_sample": {
                "data": ["0", "1", "2", "3", "4", "5"],
                "dims": ("_sample",),
                "attrs": {},
            },
        },
        "attrs": {},
        "dims": ("_sample",),
        "data": [
            -3.7523451473100002,
            -0.421128991488,
            -0.536290093143,
            -0.9076827328799999,
            -0.251889685747,
            -0.602998035829,
        ],
        "name": "phenotype",
    }

    y = DataArray.from_dict(data)

    data = {
        "coords": {
            "fid": {
                "data": [
                    "HG00111",
                    "HG00112",
                    "HG00116",
                    "HG00121",
                    "HG00133",
                    "HG00135",
                    "HG00142",
                ],
                "dims": ("sample",),
                "attrs": {},
            },
            "iid": {
                "data": [
                    "HG00111",
                    "HG00112",
                    "HG00116",
                    "HG00121",
                    "HG00133",
                    "HG00135",
                    "HG00142",
                ],
                "dims": ("sample",),
                "attrs": {},
            },
            "father": {
                "data": ["0", "0", "0", "0", "0", "0", "0"],
                "dims": ("sample",),
                "attrs": {},
            },
            "mother": {
                "data": ["0", "0", "0", "0", "0", "0", "0"],
                "dims": ("sample",),
                "attrs": {},
            },
            "gender": {
                "data": ["0", "0", "0", "0", "0", "0", "0"],
                "dims": ("sample",),
                "attrs": {},
            },
            "trait": {
                "data": ["-9", "-9", "-9", "-9", "-9", "-9", "-9"],
                "dims": ("sample",),
                "attrs": {},
            },
            "i": {"data": [0, 1], "dims": ("candidate",), "attrs": {}},
            "sample": {
                "data": [
                    "HG00111",
                    "HG00112",
                    "HG00116",
                    "HG00121",
                    "HG00133",
                    "HG00135",
                    "HG00142",
                ],
                "dims": ("sample",),
                "attrs": {},
            },
            "chrom": {"data": ["22", "22"], "dims": ("candidate",), "attrs": {}},
            "snp": {
                "data": ["rs146752890", "rs62224610"],
                "dims": ("candidate",),
                "attrs": {},
            },
            "cm": {"data": [0.0, 0.0], "dims": ("candidate",), "attrs": {}},
            "pos": {"data": [16050612, 16051347], "dims": ("candidate",), "attrs": {}},
            "a0": {"data": ["G", "C"], "dims": ("candidate",), "attrs": {}},
            "a1": {"data": ["C", "G"], "dims": ("candidate",), "attrs": {}},
            "candidate": {
                "data": ["rs146752890", "rs62224610"],
                "dims": ("candidate",),
                "attrs": {},
            },
        },
        "attrs": {},
        "dims": ("sample", "candidate"),
        "data": [
            [2.0, 0.0],
            [1.0, 2.0],
            [2.0, 2.0],
            [1.0, 2.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [2.0, 2.0],
        ],
        "name": "genotype",
    }

    G = DataArray.from_dict(data)

    data = conform_dataset(y, G=G)
    assert_equal(
        data["y"].coords["sample"][:3].values, ["HG00111", "HG00112", "HG00116"]
    )
    assert_equal(data["y"].shape, (6, 1))
    assert_equal(data["y"].dims, ("sample", "trait"))

    data = {
        "coords": {
            "trait": {"data": "gene1", "dims": (), "attrs": {}},
            "sample": {
                "data": ["0", "1", "2", "3", "4", "5"],
                "dims": ("sample",),
                "attrs": {},
            },
        },
        "attrs": {},
        "dims": ("sample",),
        "data": [
            -3.7523451473100002,
            -0.421128991488,
            -0.536290093143,
            -0.9076827328799999,
            -0.251889685747,
            -0.602998035829,
        ],
        "name": "phenotype",
    }

    y = DataArray.from_dict(data)
    data = conform_dataset(y, G=G)
    assert_equal(data["y"].shape, (0, 1))
    assert_equal(data["y"].dims, ("sample", "trait"))

    data = {
        "coords": {"trait": {"data": "gene1", "dims": (), "attrs": {}}},
        "attrs": {},
        "dims": ("sample",),
        "data": [
            -3.7523451473100002,
            -0.421128991488,
            -0.536290093143,
            -0.9076827328799999,
            -0.251889685747,
            -0.602998035829,
        ],
        "name": "phenotype",
    }
    y = DataArray.from_dict(data)
    data = conform_dataset(y, G=G)
    assert_equal(
        data["y"].coords["sample"][:3].values, ["HG00111", "HG00112", "HG00116"]
    )
    assert_equal(data["y"].shape, (6, 1))
    assert_equal(data["y"].dims, ("sample", "trait"))
