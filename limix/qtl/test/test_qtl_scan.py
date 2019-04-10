import pytest
import scipy.stats as st
from limix.qc import normalise_covariance
from limix.qtl import scan
from limix.stats import linear_kinship
from numpy import (
    argmin,
    array,
    asarray,
    concatenate,
    dot,
    eye,
    kron,
    nan,
    ones,
    reshape,
    sqrt,
    zeros,
)
from numpy.random import RandomState
from numpy.testing import assert_allclose
from pandas import DataFrame


def test_qtl_scan_three_hypotheses():
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

    mvn = st.multivariate_normal
    m = kron(A, M) @ beta + kron(A01, G) @ alpha
    Y = unvec(mvn(m, kron(C0, K) + kron(C1, eye(n))).rvs(), (n, -1))

    idx = [[0, 1], 2, [3]]
    r = scan(G, Y, idx=idx, K=K, M=M, A=A, A0=A0, A1=A1, verbose=False)
    print(r)


def test_qtl_scan_two_hypotheses():
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

    mvn = st.multivariate_normal
    m = kron(A, M) @ beta + kron(A01, G) @ alpha
    Y = unvec(mvn(m, kron(C0, K) + kron(C1, eye(n))).rvs(), (n, -1))

    idx = [[0, 1], 2, [3]]
    r = scan(G, Y, idx=idx, K=K, M=M, A=A, A1=A1, verbose=False)
    print(r)


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
    assert_allclose(pv[ix_best_snp], 1.0)
    assert_allclose(pv.values[0], 0.6684700834450028, rtol=1e-6)


def test_qtl_scan_glmm_binomial():
    random = RandomState(0)
    nsamples = 50

    X = random.randn(50, 2)
    G = random.randn(50, 100)
    K = dot(G, G.T)
    ntrials = random.randint(1, 100, nsamples)
    z = dot(G, random.randn(100)) / sqrt(100)

    successes = zeros(len(ntrials), int)
    for i, nt in enumerate(ntrials):
        for _ in range(nt):
            successes[i] += int(z[i] + 0.5 * random.randn() > 0)

    result = scan(X, successes, ("binomial", ntrials), K, verbose=False)
    pv = result.stats["pv20"]
    assert_allclose(pv, [0.4247691183134664, 0.707939424330281], atol=1e-6, rtol=1e-6)


def test_qtl_scan_glmm_wrong_dimensions():
    random = RandomState(0)
    nsamples = 50

    X = random.randn(50, 2)
    G = random.randn(50, 100)
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
    nsamples = 50

    X = random.randn(50, 2)
    G = random.randn(50, 100)
    K = dot(G, G.T)
    ntrials = random.randint(1, 2, nsamples)
    z = dot(G, random.randn(100)) / sqrt(100)

    successes = zeros(len(ntrials), int)
    for i, nt in enumerate(ntrials):
        for _ in range(nt):
            successes[i] += int(z[i] + 0.5 * random.randn() > 0)

    result = scan(X, successes, "bernoulli", K, verbose=False)
    pv = result.stats["pv20"]
    assert_allclose(pv, [0.34545196133488565, 0.3553914023889331], rtol=1e-5)


def test_qtl_scan_glmm_bernoulli_nokinship():
    random = RandomState(0)
    nsamples = 50

    X = random.randn(50, 2)
    G = random.randn(50, 100)
    ntrials = random.randint(1, 2, nsamples)
    z = dot(G, random.randn(100)) / sqrt(100)

    successes = zeros(len(ntrials), int)
    for i, nt in enumerate(ntrials):
        for _ in range(nt):
            successes[i] += int(z[i] + 0.5 * random.randn() > 0)

    result = scan(X, successes, "bernoulli", verbose=False)
    pv = result.stats["pv20"]
    assert_allclose(pv, [0.9255237731412944, 0.17024464676060735], rtol=1e-5)


def test_qtl_scan_lm():
    random = RandomState(0)
    nsamples = 50

    G = random.randn(50, 100)

    y = dot(G, random.randn(100)) / sqrt(100) + 0.2 * random.randn(nsamples)

    M = G[:, :5]
    X = G[:, 5:]
    result = scan(X, y, "normal", M=M, verbose=False)
    pv = result.stats["pv20"]
    assert_allclose(pv[:2], [0.13302121289855245, 0.031231550764835345], rtol=1e-5)


def test_qtl_scan_gmm_binomial():
    random = RandomState(0)
    nsamples = 50

    X = random.randn(nsamples, 2)
    ntrials = random.randint(1, nsamples, nsamples)
    z = dot(X, random.randn(2))

    successes = zeros(len(ntrials), int)
    for i in range(len(ntrials)):
        for _ in range(ntrials[i]):
            successes[i] += int(z[i] + 0.5 * random.randn() > 0)

    result = scan(X, successes, ("binomial", ntrials), verbose=False)
    pv = result.stats["pv20"]
    assert_allclose(pv, [0.38010406141, 1.6679817550561325e-21], rtol=1e-5, atol=1e-5)


def test_qtl_finite():
    random = RandomState(0)
    nsamples = 50

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
