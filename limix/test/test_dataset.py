from numpy import array, dtype
from numpy.random import RandomState
from numpy.testing import assert_, assert_array_equal, assert_equal
from pandas import DataFrame, Series
from xarray import DataArray

from limix._data import conform_dataset


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
    from limix._data import asarray

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

    for xi in x:
        y = asarray(xi, "trait", ["sample", "trait"])
        assert_equal(y.dtype, dtype("float64"))
        assert_array_equal(y.shape, (3, 1))
        assert_(isinstance(y, DataArray))
        if isinstance(xi, Series):
            assert_array_equal(list(xi.index), list(y.coords["sample"].values))
        if isinstance(xi, DataFrame):
            assert_array_equal(list(xi.columns), list(y.coords["trait"].values))

        is_dask = (
            hasattr(xi, "chunks")
            and xi.chunks is not None
            or hasattr(xi, "values")
            and hasattr(xi, "values")
            and hasattr(xi.values, "chunks")
            and xi.values.chunks is not None
        )

        assert_equal(is_dask, y.chunks is not None)
        assert_array_equal(np.asarray(xi).ravel(), np.asarray(y).ravel())


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
