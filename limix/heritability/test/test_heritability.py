from numpy import dot, exp, sqrt, zeros
from numpy.random import RandomState
from numpy.testing import assert_allclose
from numpy_sugar.linalg import economic_qs

from limix.heritability import estimate


def test_heritability_estimate_binomial():
    random = RandomState(0)
    nsamples = 50

    X = random.randn(50, 2)
    G = random.randn(50, 100)
    K = dot(G, G.T)
    ntrials = random.randint(1, 100, nsamples)
    z = dot(G, random.randn(100)) / sqrt(100)

    successes = zeros(len(ntrials), int)
    for i in range(len(ntrials)):
        for j in range(ntrials[i]):
            successes[i] += int(z[i] + 0.5 * random.randn() > 0)

    y = (successes, ntrials)

    assert_allclose(
        estimate(y, "binomial", K, verbose=False), 0.8937449693087052, rtol=1e-5
    )


def test_heritability_estimate_poisson():
    random = RandomState(0)
    nsamples = 50

    G = random.randn(50, 100)
    K = dot(G, G.T)
    z = dot(G, random.randn(100)) / sqrt(100)
    y = random.poisson(exp(z))

    assert_allclose(
        estimate(y, "poisson", K, verbose=False), 0.6998737749318877, rtol=1e-5
    )
