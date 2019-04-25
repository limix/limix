import pytest
from xarray import DataArray
import scipy.stats as st
from numpy import (
    argmin,
    array,
    concatenate,
    dot,
    exp,
    eye,
    kron,
    nan,
    reshape,
    sqrt,
    zeros,
)
from numpy.random import RandomState
from numpy.testing import assert_allclose, assert_array_equal
from pandas import DataFrame

from limix.qc import normalise_covariance
from limix.qtl import scan
from limix.stats import linear_kinship, multivariate_normal as mvn


def _test_qtl_scan_st(lik):
    random = RandomState(0)
    n = 30
    ncovariates = 3

    M = random.randn(n, ncovariates)

    v0 = random.rand()
    v1 = random.rand()

    G = random.randn(n, 4)

    K = random.randn(n, n + 1)
    K = normalise_covariance(K @ K.T)

    beta = random.randn(ncovariates)
    alpha = random.randn(G.shape[1])

    m = M @ beta + G @ alpha
    y = mvn(random, m, v0 * K + v1 * eye(n))

    idx = [[0, 1], 2, [3]]

    if lik == "poisson":
        y = random.poisson(exp(y))
    elif lik == "bernoulli":
        y = random.binomial(1, 1 / (1 + exp(-y)))
    elif lik == "probit":
        y = random.binomial(1, st.norm.cdf(y))
    elif lik == "binomial":
        ntrials = random.randint(0, 30, len(y))
        y = random.binomial(ntrials, 1 / (1 + exp(-y)))
        lik = (lik, ntrials)

    r = scan(G, y, lik=lik, idx=idx, K=K, M=M, verbose=False)
    str(r)
    str(r.stats.head())
    str(r.effsizes["h2"].head())
    str(r.h0.trait)
    str(r.h0.likelihood)
    str(r.h0.lml)
    str(r.h0.effsizes)
    str(r.h0.variances)


def test_qtl_scan_st():
    _test_qtl_scan_st("normal")
    _test_qtl_scan_st("poisson")
    _test_qtl_scan_st("bernoulli")
    _test_qtl_scan_st("probit")
    _test_qtl_scan_st("binomial")


def test_qtl_scan_three_hypotheses_mt():
    random = RandomState(0)
    n = 30
    ntraits = 2
    ncovariates = 3

    A = random.randn(ntraits, ntraits)
    A = A @ A.T
    M = random.randn(n, ncovariates)

    C0 = random.randn(ntraits, ntraits)
    C0 = C0 @ C0.T

    C1 = random.randn(ntraits, ntraits)
    C1 = C1 @ C1.T

    G = random.randn(n, 4)

    A0 = random.randn(ntraits, 1)
    A1 = random.randn(ntraits, 2)
    A01 = concatenate((A0, A1), axis=1)

    K = random.randn(n, n + 1)
    K = normalise_covariance(K @ K.T)

    beta = vec(random.randn(ntraits, ncovariates))
    alpha = vec(random.randn(A01.shape[1], G.shape[1]))

    m = kron(A, M) @ beta + kron(A01, G) @ alpha
    Y = unvec(mvn(random, m, kron(C0, K) + kron(C1, eye(n))), (n, -1))

    idx = [[0, 1], 2, [3]]
    r = scan(G, Y, idx=idx, K=K, M=M, A=A, A0=A0, A1=A1, verbose=False)
    str(r)


def test_qtl_scan_two_hypotheses_mt():
    random = RandomState(0)
    n = 30
    ntraits = 2
    ncovariates = 3

    A = random.randn(ntraits, ntraits)
    A = A @ A.T
    M = random.randn(n, ncovariates)

    C0 = random.randn(ntraits, ntraits)
    C0 = C0 @ C0.T

    C1 = random.randn(ntraits, ntraits)
    C1 = C1 @ C1.T

    G = random.randn(n, 4)

    A0 = random.randn(ntraits, 1)
    A1 = random.randn(ntraits, 2)
    A01 = concatenate((A0, A1), axis=1)

    K = random.randn(n, n + 1)
    K = normalise_covariance(K @ K.T)

    beta = vec(random.randn(ntraits, ncovariates))
    alpha = vec(random.randn(A01.shape[1], G.shape[1]))

    m = kron(A, M) @ beta + kron(A01, G) @ alpha
    Y = unvec(mvn(random, m, kron(C0, K) + kron(C1, eye(n))), (n, -1))

    idx = [[0, 1], 2, [3]]
    r = scan(G, Y, idx=idx, K=K, M=M, A=A, A1=A1, verbose=False)
    str(r)


def test_qtl_scan_two_hypotheses_mt_A0A1_none():
    random = RandomState(0)
    n = 30
    ntraits = 2
    ncovariates = 3

    A = random.randn(ntraits, ntraits)
    A = A @ A.T
    M = random.randn(n, ncovariates)

    C0 = random.randn(ntraits, ntraits)
    C0 = C0 @ C0.T

    C1 = random.randn(ntraits, ntraits)
    C1 = C1 @ C1.T

    G = random.randn(n, 4)

    A1 = eye(ntraits)

    K = random.randn(n, n + 1)
    K = normalise_covariance(K @ K.T)

    beta = vec(random.randn(ntraits, ncovariates))
    alpha = vec(random.randn(A1.shape[1], G.shape[1]))

    m = kron(A, M) @ beta + kron(A1, G) @ alpha
    Y = unvec(mvn(random, m, kron(C0, K) + kron(C1, eye(n))), (n, -1))
    Y = DataArray(Y, dims=["sample", "trait"], coords={"trait": ["WA", "Cx"]})

    idx = [[0, 1], 2, [3]]
    r = scan(G, Y, idx=idx, K=K, M=M, A=A, verbose=False)
    df = r.effsizes["h2"]
    df = df[df["test"] == 0]
    assert_array_equal(df["trait"], ["WA"] * 3 + ["Cx"] * 3 + [None] * 4)
    assert_array_equal(
        df["env"], [None] * 6 + ["env1_WA", "env1_WA", "env1_Cx", "env1_Cx"]
    )
    str(r)


def test_qtl_scan_lmm():
    random = RandomState(0)
    nsamples = 50

    G = random.randn(50, 100)
    K = linear_kinship(G[:, 0:80], verbose=False)

    y = dot(G, random.randn(100)) / sqrt(100) + 0.2 * random.randn(nsamples)

    M = G[:, :5]
    X = G[:, 68:70]

    result = scan(X, y, lik="normal", K=K, M=M, verbose=False)

    pv = result.stats["pv20"]

    ix_best_snp = argmin(array(pv))
    M = concatenate((M, X[:, [ix_best_snp]]), axis=1)
    result = scan(X, y, "normal", K, M=M, verbose=False)
    pv = result.stats["pv20"]
    assert_allclose(pv[ix_best_snp], 1.0, atol=1e-6)


def test_qtl_scan_lmm_nokinship():
    random = RandomState(0)
    nsamples = 50

    G = random.randn(50, 100)
    K = linear_kinship(G[:, 0:80], verbose=False)

    y = dot(G, random.randn(100)) / sqrt(100) + 0.2 * random.randn(nsamples)

    M = G[:, :5]
    X = G[:, 68:70]

    result = scan(X, y, "normal", K, M=M, verbose=False)
    pv = result.stats["pv20"].values
    assert_allclose(pv[:2], [8.159539103135342e-05, 0.10807353641893498], atol=1e-5)


def test_qtl_scan_lmm_repeat_samples_by_index():
    random = RandomState(0)
    nsamples = 30
    samples = ["sample{}".format(i) for i in range(nsamples)]

    G = random.randn(nsamples, 100)
    G = DataFrame(data=G, index=samples)

    K = linear_kinship(G.values[:, 0:80], verbose=False)
    K = DataFrame(data=K, index=samples, columns=samples)

    y0 = dot(G, random.randn(100)) / sqrt(100) + 0.2 * random.randn(nsamples)
    y1 = dot(G, random.randn(100)) / sqrt(100) + 0.2 * random.randn(nsamples)
    y = concatenate((y0, y1))
    y = DataFrame(data=y, index=samples + samples)

    M = G.values[:, :5]
    X = G.values[:, 68:70]
    M = DataFrame(data=M, index=samples)
    X = DataFrame(data=X, index=samples)

    result = scan(X, y, "normal", K, M=M, verbose=False)
    pv = result.stats["pv20"]
    assert_allclose(pv.values[0], 0.9920306566395604, rtol=1e-6)

    ix_best_snp = argmin(array(result.stats["pv20"]))

    M = concatenate((M, X.loc[:, [ix_best_snp]]), axis=1)
    M = DataFrame(data=M, index=samples)

    result = scan(X, y, "normal", K, M=M, verbose=False)
    pv = result.stats["pv20"]
    assert_allclose(pv[ix_best_snp], 1.0, rtol=1e-6)
    assert_allclose(pv.values[0], 0.6684700834450028, rtol=1e-6)


def test_qtl_scan_glmm_binomial():
    random = RandomState(0)
    nsamples = 25

    X = random.randn(nsamples, 2)
    G = random.randn(nsamples, 100)
    K = dot(G, G.T)
    ntrials = random.randint(1, 100, nsamples)
    z = dot(G, random.randn(100)) / sqrt(100)

    successes = zeros(len(ntrials), int)
    for i, nt in enumerate(ntrials):
        for _ in range(nt):
            successes[i] += int(z[i] + 0.5 * random.randn() > 0)

    result = scan(X, successes, ("binomial", ntrials), K, verbose=False)
    pv = result.stats["pv20"]
    assert_allclose(pv, [0.9315770010211236, 0.8457015828837173], atol=1e-6, rtol=1e-6)


def test_qtl_scan_glmm_wrong_dimensions():
    random = RandomState(0)
    nsamples = 25

    X = random.randn(nsamples, 2)
    G = random.randn(nsamples, 100)
    K = dot(G, G.T)
    ntrials = random.randint(1, 100, nsamples)
    z = dot(G, random.randn(100)) / sqrt(100)

    successes = zeros(len(ntrials), int)
    for i, nt in enumerate(ntrials):
        for _ in range(nt):
            successes[i] += int(z[i] + 0.5 * random.randn() > 0)

    M = random.randn(49, 2)
    scan(X, successes, ("binomial", ntrials), K, M=M, verbose=False)


def test_qtl_scan_glmm_bernoulli():
    random = RandomState(0)
    nsamples = 25

    X = random.randn(nsamples, 2)
    G = random.randn(nsamples, 100)
    K = dot(G, G.T)
    ntrials = random.randint(1, 2, nsamples)
    z = dot(G, random.randn(100)) / sqrt(100)

    successes = zeros(len(ntrials), int)
    for i, nt in enumerate(ntrials):
        for _ in range(nt):
            successes[i] += int(z[i] + 0.5 * random.randn() > 0)

    result = scan(X, successes, "bernoulli", K, verbose=False)
    pv = result.stats["pv20"]
    assert_allclose(pv, [0.3399326545917558, 0.8269454251659921], rtol=1e-5)


def test_qtl_scan_glmm_bernoulli_nokinship():
    random = RandomState(0)
    nsamples = 25

    X = random.randn(nsamples, 2)
    G = random.randn(nsamples, 100)
    ntrials = random.randint(1, 2, nsamples)
    z = dot(G, random.randn(100)) / sqrt(100)

    successes = zeros(len(ntrials), int)
    for i, nt in enumerate(ntrials):
        for _ in range(nt):
            successes[i] += int(z[i] + 0.5 * random.randn() > 0)

    result = scan(X, successes, "bernoulli", verbose=False)
    pv = result.stats["pv20"]
    assert_allclose(pv, [0.3399067917883736, 0.8269568797830423], rtol=1e-5)


def test_qtl_scan_lm():
    random = RandomState(0)
    nsamples = 25

    G = random.randn(nsamples, 100)

    y = dot(G, random.randn(100)) / sqrt(100) + 0.2 * random.randn(nsamples)

    M = G[:, :5]
    X = G[:, 5:]
    result = scan(X, y, "normal", M=M, verbose=False)
    pv = result.stats["pv20"]
    assert_allclose(pv[:2], [0.02625506841465465, 0.9162689001409643], rtol=1e-5)


def test_qtl_scan_gmm_binomial():
    random = RandomState(0)
    nsamples = 25

    X = random.randn(nsamples, 2)
    ntrials = random.randint(1, nsamples, nsamples)
    z = dot(X, random.randn(2))

    successes = zeros(len(ntrials), int)
    for i in range(len(ntrials)):
        for _ in range(ntrials[i]):
            successes[i] += int(z[i] + 0.5 * random.randn() > 0)

    result = scan(X, successes, ("binomial", ntrials), verbose=False)
    pv = result.stats["pv20"]
    assert_allclose(
        pv, [2.4604711379400065e-06, 0.01823278752006871], rtol=1e-5, atol=1e-5
    )


def test_qtl_finite():
    random = RandomState(0)
    nsamples = 20

    X = random.randn(50, 2)
    G = random.randn(50, 100)
    K = dot(G, G.T)
    ntrials = random.randint(1, 100, nsamples)
    z = dot(G, random.randn(100)) / sqrt(100)

    successes = zeros(len(ntrials), int)
    for i, nt in enumerate(ntrials):
        for _ in range(nt):
            successes[i] += int(z[i] + 0.5 * random.randn() > 0)

    successes = successes.astype(float)
    ntrials = ntrials.astype(float)

    successes[0] = nan
    with pytest.raises(ValueError):
        scan(X, successes, ("binomial", ntrials), K, verbose=False)
    successes[0] = 1.0

    K[0, 0] = nan
    with pytest.raises(ValueError):
        scan(X, successes, ("binomial", ntrials), K, verbose=False)
    K[0, 0] = 1.0

    X[0, 0] = nan
    with pytest.raises(ValueError):
        scan(X, successes, ("binomial", ntrials), K, verbose=False)
    X[0, 0] = 1.0


def vec(x):
    return reshape(x, (-1,) + x.shape[2:], order="F")


def unvec(x, shape):
    return reshape(x, shape, order="F")
