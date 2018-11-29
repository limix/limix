from numpy import dot, sqrt
from numpy.random import RandomState
from numpy.testing import assert_allclose

from limix.qtl import st_scan


def test_qtl_lm():
    random = RandomState(0)
    nsamples = 50

    G = random.randn(50, 100)

    y = dot(G, random.randn(100)) / sqrt(100) + 0.2 * random.randn(nsamples)

    M = G[:, :5]
    X = G[:, 5:]
    model = st_scan(X, y, "normal", M=M, verbose=False)
    pv = model.variant_pvalues
    assert_allclose(pv[:2], [0.133021212899, 0.0312315507648], rtol=1e-4)
