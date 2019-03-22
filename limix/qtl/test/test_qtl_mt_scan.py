from limix.qtl import mt_scan
from numpy.random import RandomState
from numpy.testing import assert_allclose


def test_qtl_mt_scan():
    random = RandomState(0)
    n = 5
    ntraits = 3
    Y = random.randn(n, ntraits)
    G = random.randn(n, 4)
    K = random.randn(n, n)
    K = K @ K.T
    K /= K.diagonal().mean()

    r = mt_scan(G, Y, K=K, verbose=False)
    print(r.null_effsizes())

    # A = random.randn(3, 3)
    # A = A @ A.T
    # M = random.randn(n, 2)
