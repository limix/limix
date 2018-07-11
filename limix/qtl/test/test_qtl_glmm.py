from numpy import dot, sqrt, zeros
from numpy.random import RandomState
from numpy.testing import assert_allclose
from numpy_sugar.linalg import economic_qs

from limix.qtl import qtl_test_glmm


def test_qtl_glmm_binomial():
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

    lmm = qtl_test_glmm(X, y, "binomial", K, verbose=False)
    pv = lmm.getPv()
    assert_allclose(pv, [[0.44244291380638878, 0.67987701835252734]], rtol=1e-5)
