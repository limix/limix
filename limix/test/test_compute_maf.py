from limix.qc import compute_maf
from numpy.random import RandomState
from numpy.testing import assert_allclose


def test_compute_maf():
    random = RandomState(0)

    X = random.randint(0, 3, size=(100, 10))
    assert_allclose(
        compute_maf(X), [0.49, 0.49, 0.445, 0.495, 0.5, 0.45, 0.48, 0.48, 0.47, 0.435]
    )
