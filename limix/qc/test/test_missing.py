import dask.array as da
from numpy import nan
from numpy.random import RandomState
from numpy.testing import assert_array_equal

from limix.qc import count_missingness


def test_count_missingness():
    random = RandomState(0)
    X = random.randn(3, 4)
    X[0, 0] = nan
    X[0, 2] = nan
    X[2, 2] = nan
    X[:, 3] = nan
    assert_array_equal(count_missingness(X), [1, 0, 2, 3])
    X = da.from_array(X, 2)
    assert_array_equal(count_missingness(X), [1, 0, 2, 3])
