from numpy import dot, sqrt, asarray, ones, zeros
from numpy.random import RandomState
from numpy.testing import assert_allclose, assert_array_equal

from limix.qtl import st_scan, st_iscan
from limix.stats import linear_kinship


def test_qtl_st_sscan():

    random = RandomState(0)
    nsamples = 500

    G = random.randn(nsamples, 100)
    K = linear_kinship(G[:, 5:80], verbose=False)

    beta = asarray([-1, 1, 1], float)
    y = dot(G[:, :3], beta) / sqrt(3) + 5.0 * random.randn(nsamples)

    M = G[:, 10:12]
    X = G[:, :3]

    r0 = st_scan(X, y, "normal", K, M=M, verbose=False)
    E0 = zeros((nsamples, 1))
    E1 = ones((nsamples, 1))
    r1 = st_iscan(X, y, M=M, E1=E1, E0=E0, K=K, verbose=False)

    assert_allclose(r1.stats["null lml"], r0.stats["null lml"])
    assert_allclose(r1.stats["alt lml"], r0.stats["alt lml"])
    assert_allclose(r1.foreground_variance, r0.foreground_variance)
    assert_allclose(r1.background_variance, r0.background_variance)
    assert_allclose(r1.covariate_effsizes, r0.covariate_effsizes)

    effsizes = r1.alt_effsizes[r1.alt_effsizes["inter"] == "inter1_0"]["effsize"]
    assert_allclose(effsizes, r0.alt_effsizes["effsize"])

    effsizes = r1.alt_effsizes[r1.alt_effsizes["inter"] == "inter1_0"]["effsize se"]
    assert_allclose(effsizes, r0.alt_effsizes["effsize se"])

    assert_array_equal(r1.alt_effsizes["test"], [0, 0, 1, 1, 2, 2])
    assert_array_equal(r1.alt_effsizes["candidate"], ["0", "0", "1", "1", "2", "2"])

    assert_array_equal(r1.null_effsizes["test"], [0, 1, 2])
    assert_array_equal(r1.null_effsizes["candidate"], ["0", "1", "2"])
    assert_array_equal(r1.null_effsizes["effsize"], [0.0] * 3)
    assert_array_equal(r1.null_effsizes["effsize se"], [0.0] * 3)

    assert_allclose(r1.background_variance, r0.background_variance)
    assert_allclose(r1.foreground_variance, r0.foreground_variance)

    E0 = random.randn(nsamples, 2)
    E1 = random.randn(nsamples, 3)

    r = st_iscan(X, y, idx=[slice(0, 2), 2], M=M, E1=E1, E0=E0, K=K, verbose=False)
    assert_allclose(
        r.covariate_effsizes["effsizes"], [0.5993381138115996, -0.13135620851143137]
    )
    assert_allclose(r.null_effsizes["test"], [0, 0, 0, 0, 1, 1])
    assert_array_equal(r.null_effsizes["candidate"], ["0", "0", "1", "1", "2", "2"])
    assert_array_equal(
        r.null_effsizes["inter"],
        ["inter0_0", "inter0_1", "inter0_0", "inter0_1", "inter0_0", "inter0_1"],
    )
    assert_allclose(
        r.null_effsizes["effsize"],
        [
            0.18808905570133957,
            0.10414389776928647,
            -0.23532951492536575,
            -0.09804663988905885,
            0.10902383683518363,
            0.1389349173733463,
        ],
    )
    assert_allclose(
        r.null_effsizes["effsize se"],
        [
            0.1378142439333301,
            0.07630700509302424,
            0.17242767822777486,
            0.07183949909328986,
            0.140026505737559,
            0.17844327964843626,
        ],
    )
    assert_array_equal(
        r.alt_effsizes["inter"],
        ["inter0_0", "inter0_1", "inter1_0", "inter1_1", "inter1_2"] * 3,
    )

    r = st_iscan(X, y, idx=[slice(0, 2), 2], M=M, E1=E1, K=K, verbose=False)
    assert_array_equal(r.null_effsizes["inter"], ["inter0_0"] * 3)
    assert_allclose(
        r.null_effsizes["effsize"],
        [-0.7789846307039473, 0.36202092757221177, 0.47801355168983917],
    )


if __name__ == "__main__":
    test_qtl_st_sscan()
