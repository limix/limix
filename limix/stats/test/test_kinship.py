import dask.array as da
from numpy.random import RandomState
from numpy.testing import assert_allclose

from limix.stats import linear_kinship


def test_kinship_estimation():
    random = RandomState(0)
    X = random.randn(30, 40)

    K0 = linear_kinship(X, verbose=False)

    X = da.from_array(X, chunks=(5, 13))
    K1 = linear_kinship(X, verbose=False)
    assert_allclose(K0, K1)
