from __future__ import division

from dask.array import from_array
from numpy import logical_not
from numpy.random import RandomState
from numpy.testing import assert_allclose, assert_equal

from limix.stats import indep_pairwise, maf


def test_preprocess_indep_pairwise():
    random = RandomState(0)

    X = random.randn(3, 100)

    head = [True, True, False, True, False]
    tail = [True, True, False, False]

    assert_equal(indep_pairwise(X, 4, 2, 0.5, verbose=False)[:5], head)
    assert_equal(indep_pairwise(X, 4, 2, 0.5, verbose=False)[-4:], tail)

    X = from_array(X, chunks=(2, 10))
    assert_equal(indep_pairwise(X, 4, 2, 0.5, verbose=False)[:5], head)
    assert_equal(indep_pairwise(X, 4, 2, 0.5, verbose=False)[-4:], tail)


def test_preprocess_maf():
    random = RandomState(0)

    X = random.randint(0, 3, size=(100, 10))
    assert_allclose(
        maf(X), [0.49, 0.49, 0.445, 0.495, 0.5, 0.45, 0.48, 0.48, 0.47, 0.435])
