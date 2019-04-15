from numpy import dot, exp, sqrt, zeros
from numpy.random import RandomState
from numpy.testing import assert_allclose

from limix.her import estimate


def test_heritability_estimate_binomial():
    random = RandomState(0)
    nsamples = 50

    G = random.randn(nsamples, 100)
    K = dot(G, G.T)
    ntrials = random.randint(1, 100, nsamples)
    z = dot(G, random.randn(100)) / sqrt(100)

    successes = zeros(len(ntrials), int)
    for i, nt in enumerate(ntrials):
        for _ in range(nt):
            successes[i] += int(z[i] + 0.5 * random.randn() > 0)

    assert_allclose(
        estimate(successes, ("binomial", ntrials), K, verbose=False),
        0.9999992083233082,
        rtol=1e-3,
    )


def test_heritability_estimate_poisson():
    random = RandomState(0)

    G = random.randn(50, 100)
    K = dot(G, G.T)
    z = dot(G, random.randn(100)) / sqrt(100)
    y = random.poisson(exp(z))

    assert_allclose(
        estimate(y, "poisson", K, verbose=False), 0.991766763337491, rtol=1e-3
    )


def test_heritability_estimate_normal():
    random = RandomState(0)

    G = random.randn(50, 100)
    K = dot(G, G.T)
    z = dot(G, random.randn(100)) / sqrt(100)
    y = z + 0.2 * random.randn(50)

    h2 = estimate(y, "normal", K, verbose=False)
    assert_allclose(h2, 0.9751053293095242, rtol=1e-5)
