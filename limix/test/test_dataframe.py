import pytest
from numpy import array
from numpy.random import RandomState
from numpy.testing import assert_, assert_array_equal
from pandas import DataFrame, Index

from limix.dataframe import (
    _infer_samples_index,
    assert_compatible_samples_arrays,
    normalise_dataset,
)


def test_dataframe_arrays_compatibility():
    y = array([-1.2, 3.4, 0.1])
    assert_(assert_compatible_samples_arrays([]))
    assert_(assert_compatible_samples_arrays(y))
    assert_(assert_compatible_samples_arrays([y]))

    random = RandomState(0)

    K = random.randn(3, 4)
    K = K.dot(K.T)
    assert_(assert_compatible_samples_arrays([y, K]))

    M = random.randn(3, 2)
    assert_(assert_compatible_samples_arrays([y, K, M]))
    Me = random.randn(2, 2)
    assert_(not assert_compatible_samples_arrays([y, K, Me]))


def test_dataframe_samples_index_infer():
    y = array([-1.2, 3.4, 0.1])

    random = RandomState(0)

    K = random.randn(3, 4)
    K = K.dot(K.T)

    M = random.randn(3, 2)
    Me = random.randn(2, 2)

    with pytest.raises(ValueError):
        _infer_samples_index([y, K, Me])

    candidates = Index(["candidate0", "candidate1", "candidate2"])
    assert_array_equal(_infer_samples_index(y), candidates)
    assert_array_equal(_infer_samples_index([y]), candidates)
    assert_array_equal(_infer_samples_index([y, K]), candidates)
    assert_array_equal(_infer_samples_index([y, K, M]), candidates)

    df_y = DataFrame(y, index=["c0", "c1", "c2"])
    df_K = DataFrame(K, index=["c0", "c1", "c2"])
    df_Me = DataFrame(Me, index=["c0", "c1"])
    c = Index(["c0", "c1"])
    assert_array_equal(_infer_samples_index([df_y, df_K, df_Me]), c)

    cindex = Index(["c0", "c1", "c2"])
    assert_array_equal(_infer_samples_index([df_y, df_K]), cindex)
    assert_array_equal(_infer_samples_index([df_y, df_K, M]), cindex)

    df_K = DataFrame(K, index=["c2", "c0", "c1"])
    with pytest.raises(ValueError):
        _infer_samples_index([df_y, df_K, M])


def test_dataframe_normalise_dataset():
    y = array([-1.2, 3.4, 0.1])
    samples = ["sample{}".format(i) for i in range(len(y))]
    y = DataFrame(data=y, index=samples)

    random = RandomState(0)

    K = random.randn(3, 4)
    K = K.dot(K.T)
    K = DataFrame(data=K, index=samples, columns=samples)

    M = random.randn(3, 2)
    M = DataFrame(data=M, index=samples)

    data = normalise_dataset(y, "normal", M=M, K=K)
    assert_array_equal(y.values, data["y"].values)

    y = array([-1.2, 3.4, 0.1, 0.1, 0.0, -0.2])
    y = DataFrame(data=y, index=samples + samples)
    data = normalise_dataset(y, "normal", M=M, K=K)
    y = data["y"]
    M = data["M"]
    K = data["K"]
    assert_(y.shape[0] == M.shape[0])
    assert_(y.shape[0] == K.shape[0])
    assert_(y.shape[0] == K.shape[1])
    k = K.loc["sample0", :].loc[:, "sample0"].values
    assert_array_equal(k[0, 0], k[0, 1])
    assert_array_equal(k[0, 0], k[1, 1])
    assert_array_equal(k[0, 0], k[1, 0])

    with pytest.raises(ValueError):
        data = normalise_dataset(y, "normal", M=M, K=K, G=M.values.copy())

    G = M.copy()
    data = normalise_dataset(y, "normal", M=M, K=K, G=G)
    y = data["y"]
    M = data["M"]
    K = data["K"]
    G = data["G"]
    assert_(y.shape[0] == G.shape[0])
