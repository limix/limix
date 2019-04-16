from numpy import ones, zeros
from numpy.random import RandomState
from numpy.testing import assert_allclose

from limix.qc import remove_dependent_cols


def test_remove_dependent_cols():
    random = RandomState(0)

    X = random.randn(4, 5)
    X[:, 2] = 3 * X[:, 1]
    D = X[:, [0, 1, 3, 4]]
    assert_allclose(remove_dependent_cols(X), D)

    X = random.randn(4, 1)
    R = X.copy()
    assert_allclose(remove_dependent_cols(X), R)

    X = ones((3, 4))
    R = ones((3, 1))
    assert_allclose(remove_dependent_cols(X), R)

    X = zeros((3, 4))
    R = zeros((3, 0))
    assert_allclose(remove_dependent_cols(X), R)
