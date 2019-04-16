from dask.array import from_array
from numpy.random import RandomState
from numpy.testing import assert_equal

from limix.qc import indep_pairwise


def test_indep_pairwise():
    random = RandomState(0)

    X = random.randn(3, 10)

    head = [True, True, False, False, False]
    tail = [False, True, False, False]

    assert_equal(indep_pairwise(X, 4, 2, 0.5, verbose=False)[:5], head)
    assert_equal(indep_pairwise(X, 4, 2, 0.5, verbose=False)[-4:], tail)

    X = from_array(X, chunks=(2, 10))
    assert_equal(indep_pairwise(X, 4, 2, 0.5, verbose=False)[:5], head)
    assert_equal(indep_pairwise(X, 4, 2, 0.5, verbose=False)[-4:], tail)
