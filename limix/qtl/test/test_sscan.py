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
            0.367504189667873,
            0.16787513728846237,
            0.5471313168464726,
            0.847625281563343,
            0.47885258328495495,
            0.9303397085618881,
            0.2820898680960877,
            0.4345680087053413,
            0.8336863244890175,
            0.2825636286402263,
            0.9151742153205723,
        ],
        rtol=1e-5,
    )
    assert_allclose(
        pvi,
        [
            0.23201754963764076,
            0.08834969611344345,
            0.4132267634495471,
            0.6210097240361562,
            0.4956162463895071,
            0.7673058550911451,
            0.15997052199126327,
            0.2573438296290731,
            0.6935757798342528,
            0.16321585906881486,
            0.7414986547945575,
        ],
        rtol=1e-5,
    )
