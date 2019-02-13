from dask.array import from_array
from limix.qc import mean_standardize
from numpy.random import RandomState
from numpy.testing import assert_allclose


def test_mean_standardize():
    random = RandomState(0)
    X = random.randn(5, 3)

    def assert_correct(X):
        Y = mean_standardize(X)
        assert_allclose(Y.mean(), 0, atol=1e-7)
        assert_allclose(Y.std(), 1, atol=1e-7)

        for axis in [0, 1]:
            Y = mean_standardize(X, axis=axis)
            i = (axis + 1) % 2
            assert_allclose(Y.mean(axis), [0] * X.shape[i], atol=1e-7)
            assert_allclose(Y.std(axis), [1] * X.shape[i], atol=1e-7)

    assert_correct(X)
    assert_correct(from_array(X, chunks=(2, 1)))
