import dask.array as da
from numpy.random import RandomState
from numpy.testing import assert_array_equal

from limix.qc import unique_variants


def test_unique_variants():
    random = RandomState(0)
    X = random.randn(3, 3)
    X[:, 2] = X[:, 0]
    assert_array_equal(unique_variants(X), X[:, 1:])
    X = da.from_array(X, 1)
    assert_array_equal(unique_variants(X), X[:, 1:])
