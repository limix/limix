from numpy.random import RandomState
from numpy.testing import assert_allclose

from limix.stats import Chi2Mixture


def test_chi2mixture():
    dof = 2
    mixture = 0.2
    n = 100

    random = RandomState(1)
    x = random.chisquare(dof, n)
    n0 = int((1 - mixture) * n)
    idxs = random.choice(n, n0, replace=False)
    x[idxs] = 0

    chi2mix = Chi2Mixture(
        scale_min=0.1,
        scale_max=5.0,
        dof_min=0.1,
        dof_max=5.0,
        qmax=0.1,
        tol=4e-3,
        lrt=x,
    )
    chi2mix.estimate_chi2mixture(x)
    pv = chi2mix.sf([0.0, 0.2])
    assert_allclose(pv, [0.19999999999999996, 0.1412935752078675])
    assert_allclose(chi2mix.scale, 1.9808080808080812)
    assert_allclose(chi2mix.dof, 0.891919191919192)
    assert_allclose(chi2mix.mixture, 0.199999999999999960)
