from numpy.random import RandomState
from numpy.testing import assert_allclose

from limix.stats import empirical_pvalues, multipletests


def test_empirical_pvalues():

    random = RandomState(1)

    x0 = random.chisquare(1, 1000)
    x1 = random.chisquare(1, 10000)

    pv = empirical_pvalues(x0, x1)

    assert_allclose(pv[:3], [0.5599, 1.0, 0.8389])
    assert_allclose(pv[-3:], [0.249, 0.3278, 0.4848])


def test_multipletests():
    random = RandomState(5)

    pv = multipletests(random.chisquare(1, 10))[1]
    assert_allclose(
        pv,
        [
            0.56399935,
            0.55249248,
            0.96090415,
            0.99606969,
            0.74247325,
            0.12299317,
            1.0,
            0.91970786,
            0.99606969,
            0.99595594,
        ],
    )
