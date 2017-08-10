import limix
from numpy.random import RandomState
from numpy.testing import assert_allclose
import dask.array as da


def test_qc_trans_mean_standardize():
    random = RandomState(0)
    X = random.randn(5, 3)

    def assert_correct(X):
        Y = limix.qc.mean_standardize(X)
        assert_allclose(Y.mean(), 0, atol=1e-7)
        assert_allclose(Y.std(), 1, atol=1e-7)

        for axis in [0, 1]:
            Y = limix.qc.mean_standardize(X, axis=axis)
            i = (axis + 1) % 2
            assert_allclose(Y.mean(axis), [0] * X.shape[i], atol=1e-7)
            assert_allclose(Y.std(axis), [1] * X.shape[i], atol=1e-7)

    assert_correct(X)
    assert_correct(da.from_array(X, chunks=(2, 1)))
