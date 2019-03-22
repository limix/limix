from numpy.random import RandomState
from numpy.testing import assert_array_equal, assert_allclose
from xarray import concat

from limix._data import asarray as _asarray, conform_dataset
from limix.qtl._result import ScanResultFactory


def test_qtl_result():

    random = RandomState(0)
    n = 5
    ntraits = 3
    Y = random.randn(n, ntraits)
    G = random.randn(n, 20)
    K = random.randn(n, n)
    K = K @ K.T
    K /= K.diagonal().mean()
    M = random.randn(n, 2)
    data = conform_dataset(Y, M, G=G, K=K)

    Y = data["y"]
    M = data["M"]
    G = data["G"]
    K = data["K"]

    A0 = random.randn(ntraits, 2)
    A1 = random.randn(ntraits, 3)

    A0 = _asarray(A0, "env0", ["sample", "env"])
    A1 = _asarray(A1, "env1", ["sample", "env"])
    A01 = concat([A0, A1], dim="env")

    C0 = random.randn(ntraits, ntraits)
    C0 = C0 @ C0.T
    C1 = random.randn(ntraits, ntraits)
    C1 = C1 @ C1.T

    r = ScanResultFactory(Y.trait, M.covariate, G.candidate, A0.env, A1.env)
    r.set_h0(-10.0, random.randn(ntraits, M.shape[1]), C0, C1)

    idx = [0, 1]
    covariate_effsizes = random.randn(ntraits, len(M.covariate))
    candidate_effsizes = random.randn(ntraits, len(G[:, idx].candidate) * len(A0.env))
    h1 = {
        "lml": -9.9,
        "covariate_effsizes": covariate_effsizes,
        "candidate_effsizes": candidate_effsizes,
    }
    covariate_effsizes = random.randn(ntraits, len(M.covariate))
    candidate_effsizes = random.randn(ntraits, len(G[:, idx].candidate) * len(A01.env))
    h2 = {
        "lml": -9.88,
        "covariate_effsizes": covariate_effsizes,
        "candidate_effsizes": candidate_effsizes,
    }
    r.add_test(idx, h1, h2)

    idx = [2, 3]
    covariate_effsizes = random.randn(ntraits, len(M.covariate))
    candidate_effsizes = random.randn(ntraits, len(G[:, idx].candidate) * len(A0.env))
    h1 = {
        "lml": -9.3,
        "covariate_effsizes": covariate_effsizes,
        "candidate_effsizes": candidate_effsizes,
    }
    covariate_effsizes = random.randn(ntraits, len(M.covariate))
    candidate_effsizes = random.randn(ntraits, len(G[:, idx].candidate) * len(A01.env))
    h2 = {
        "lml": -0.1,
        "covariate_effsizes": covariate_effsizes,
        "candidate_effsizes": candidate_effsizes,
    }
    r.add_test(idx, h1, h2)

    r = r.create()

    assert_array_equal(r.effsizes["h0"].iloc[2]["effect_type"], "covariate")
    assert_array_equal(r.effsizes["h0"].iloc[2]["effect_name"], "0")
    assert_array_equal(r.effsizes["h0"].iloc[4]["trait"], "2")
    assert_array_equal(
        r.effsizes["h0"].columns, ["trait", "effect_type", "effect_name", "effsize"]
    )
    assert_allclose(r.effsizes["h0"].iloc[4]["effsize"], 0.052165079260974405)

    assert_array_equal(r.effsizes["h1"].iloc[6][0], 0)
    assert_array_equal(
        r.effsizes["h1"].iloc[6][1:-1], ["0", "candidate", "0", "env0_0"]
    )
    assert_allclose(r.effsizes["h1"].iloc[6][-1], -0.1715463312222481)
    assert_array_equal(
        r.effsizes["h1"].columns,
        ["test", "trait", "effect_type", "effect_name", "env", "effsize"],
    )

    assert_array_equal(r.effsizes["h2"].iloc[23][0], 0)
    assert_array_equal(
        r.effsizes["h2"].iloc[23][1:-1], ["0", "candidate", "1", "env1_2"]
    )
    assert_allclose(r.effsizes["h2"].iloc[23][-1], -0.4635959746460942)
    assert_array_equal(
        r.effsizes["h2"].columns,
        ["test", "trait", "effect_type", "effect_name", "env", "effsize"],
    )
    assert_array_equal(r.effsizes["h2"].iloc[1][1:-2], ["0", "covariate", "1"])
    assert_array_equal(r.effsizes["h2"].iloc[1][-2], None)

    assert_array_equal(
        r.stats.columns,
        ["lml0", "lml1", "lml2", "dof10", "dof20", "dof21", "pv10", "pv20", "pv21"],
    )
    assert_allclose(
        r.stats.iloc[0],
        [
            -10.0,
            -9.9,
            -9.88,
            12.0,
            30.0,
            18.0,
            0.9999999987251013,
            0.9999999999999998,
            0.9999999999999998,
        ],
    )
