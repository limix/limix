from limix.qtl import st_sscan
from numpy.testing import assert_allclose
from numpy.random import RandomState


def test_qtl_st_sscan():

    n = 5
    random = RandomState(0)

    G = random.randn(n, 1)
    y = random.randn(n, 1)
    E = random.randn(n, 5)
    E *= E

    r = st_sscan(G, y, E, tests=["inter", "assoc"], verbose=False)
    assert_allclose(
        [r.pvi.item(), r.pva.item()], [0.40750599278681343, 0.4636304356063339]
    )
