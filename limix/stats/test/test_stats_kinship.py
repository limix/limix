import dask.array as da
from numpy.testing import assert_allclose
from numpy.random import RandomState
from limix.stats import linear_kinship


def test_stats_kinship_estimation():
    random = RandomState(0)
    X = random.randn(30, 40)

    K0 = linear_kinship(X, verbose=False)

    X = da.from_array(X, chunks=(5, 13))
    K1 = linear_kinship(X, verbose=False)
    assert_allclose(K0, K1)
