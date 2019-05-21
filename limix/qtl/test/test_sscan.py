from numpy.testing import assert_allclose


def test_qtl_sscan():
    import numpy as np

    from limix.qtl import sscan

    G = np.load("limix/qtl/test/G.npy")
    y = np.load("limix/qtl/test/y.npy")
    E = np.load("limix/qtl/test/E.npy")
    pvi, pva = sscan(G, y, E, tests=["inter", "assoc"], verbose=False)
    assert_allclose(
        pva,
        [
            0.3675160630221894,
            0.16788187986066305,
            0.5471294447654433,
            0.8476346909821835,
            0.4788508197432795,
            0.9303387337939423,
            0.2821001447449186,
            0.4345808638038694,
            0.8336859133033261,
            0.28257408571730547,
            0.9151808673835736,
        ],
        rtol=1e-5,
    )
    assert_allclose(
        pvi,
        [
            0.23200536259670557,
            0.08834969611344368,
            0.4132267634495468,
            0.6209687434734553,
            0.49561574770791705,
            0.7673032291151423,
            0.1599723078973464,
            0.2573458535281681,
            0.6935792867803385,
            0.16321675339586683,
            0.7414911089556089,
        ],
        rtol=1e-5,
    )
