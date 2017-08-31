from __future__ import division

import pytest
from numpy import dot, nan, sqrt, zeros
from numpy.random import RandomState
from numpy.testing import assert_allclose

from limix.qtl import scan


def test_qtl_finite():
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

    y = (successes.astype(float), ntrials.astype(float))

    y[0][0] = nan
    with pytest.raises(ValueError):
        scan(X, y, 'binomial', K, verbose=False)
    y[0][0] = 1.0

    K[0, 0] = nan
    with pytest.raises(ValueError):
        scan(X, y, 'binomial', K, verbose=False)
    K[0, 0] = 1.0

    X[0, 0] = nan
    with pytest.raises(ValueError):
        scan(X, y, 'binomial', K, verbose=False)
    X[0, 0] = 1.0
