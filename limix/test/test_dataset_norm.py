import pytest
from numpy import array
from numpy.random import RandomState
from numpy.testing import assert_array_equal, assert_equal
from pandas import DataFrame

from limix.dataset_norm import normalise_dataset, _infer_samples_index


def test_dataset_norm_samples_index_infer():
    y = array([-1.2, 3.4, 0.1])

    random = RandomState(0)

    K = random.randn(3, 4)
    K = K.dot(K.T)

    M = random.randn(3, 2)
    Me = random.randn(2, 2)

    # with pytest.raises(ValueError):
    #     _infer_samples_index([y, K, Me])

    candidates = [0, 1, 2]
    assert_array_equal(_infer_samples_index([y]), candidates)
    assert_array_equal(_infer_samples_index([y, K]), candidates)
    assert_array_equal(_infer_samples_index([y, K, M]), candidates)

    df_y = DataFrame(y, index=["c0", "c1", "c2"])
    df_K = DataFrame(K, index=["c0", "c1", "c2"])
    df_Me = DataFrame(Me, index=["c0", "c1"])
    c = ["c0", "c1"]

    assert_array_equal(_infer_samples_index([df_y, df_K, df_Me]), c)

    cindex = ["c0", "c1", "c2"]
    assert_array_equal(_infer_samples_index([df_y, df_K]), cindex)
    assert_array_equal(_infer_samples_index([df_y, df_K, M]), cindex)

    df_K = DataFrame(K, index=["c2", "c0", "c1"])
    with pytest.raises(ValueError):
        _infer_samples_index([df_y, df_K, M])


def test_dataset_norm_normalise_dataset():
    y = array([-1.2, 3.4, 0.1])
    samples = ["sample{}".format(i) for i in range(len(y))]
    y = DataFrame(data=y, index=samples)

    random = RandomState(0)

    K = random.randn(3, 4)
    K = K.dot(K.T)
    K = DataFrame(data=K, index=samples, columns=samples)

    M = random.randn(3, 2)
    M = DataFrame(data=M, index=samples)

    data = normalise_dataset(y, M=M, K=K)
    assert_array_equal(y.values, data["y"].values)

    y = array([-1.2, 3.4, 0.1, 0.1, 0.0, -0.2])
    with pytest.raises(ValueError):
        data = normalise_dataset(DataFrame(data=y, index=samples + samples), M=M, K=K)

    y = DataFrame(data=y[: len(samples)], index=samples)
    y = data["y"]
    M = data["M"]
    K = data["K"]

    assert_equal(y.shape[0], M.shape[0])
    assert_equal(y.shape[0], K.shape[0])
    assert_equal(y.shape[0], K.shape[1])

    G = M.copy()
    data = normalise_dataset(y, M=M, K=K, G=G)
    y = data["y"]
    M = data["M"]
    K = data["K"]
    G = data["G"]
    assert_equal(y.shape[0], G.shape[0])


def main():
    test_dataset_norm_samples_index_infer()
    test_dataset_norm_normalise_dataset()


if __name__ == "__main__":
    main()
