from numpy.random import RandomState
from numpy.testing import assert_allclose

from limix.stats import allele_expectation, allele_frequency, compute_dosage


def test_allele_frequency():
    X = RandomState(0).randn(5, 3)
    assert_allclose(
        allele_frequency(X),
        [2.0422233965933585, 0.7940255259577332, 0.5987926640480571],
    )


def test_compute_dosage():
    X = RandomState(0).randn(2, 3, 4)
    assert_allclose(
        compute_dosage(X).tolist(),
        [
            [2.240893199201458, -0.1513572082976979, 1.454273506962975],
            [0.33367432737426683, -0.8540957393017248, -0.7421650204064419],
        ],
    )


def test_allele_expectation():
    p = RandomState(0).rand(2, 3)
    p = p / p.sum(1)[:, None]
    assert_allclose(
        allele_expectation(p, 2, 2),
        [
            [0.9710998244964055, 1.028900175503594],
            [0.9374325310073951, 1.0625674689926048],
        ],
    )
