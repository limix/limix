from numpy import concatenate, dot, sqrt, argmin, array, zeros
from numpy.random import RandomState
from numpy.testing import assert_allclose
from pandas import DataFrame

from limix.qtl import st_scan
from limix.stats import linear_kinship


def test_qtl_st_scan_lmm():
    random = RandomState(0)
    nsamples = 50

    G = random.randn(50, 100)
    K = linear_kinship(G[:, 0:80], verbose=False)

    y = dot(G, random.randn(100)) / sqrt(100) + 0.2 * random.randn(nsamples)

    M = G[:, :5]
    X = G[:, 68:70]

    result = st_scan(X, y, "normal", K, M=M, verbose=False)
    pv = result.stats["pvalue"]

    ix_best_snp = argmin(array(result.stats["pvalue"]))

    M = concatenate((M, X[:, [ix_best_snp]]), axis=1)

    result = st_scan(X, y, "normal", K, M=M, verbose=False)
    pv = result.stats["pvalue"]
    assert_allclose(pv[ix_best_snp], 1.0)


def test_qtl_st_scan_lmm_nokinship():
    random = RandomState(0)
    nsamples = 50

    G = random.randn(50, 100)
    K = linear_kinship(G[:, 0:80], verbose=False)

    y = dot(G, random.randn(100)) / sqrt(100) + 0.2 * random.randn(nsamples)

    M = G[:, :5]
    X = G[:, 68:70]

    result = st_scan(X, y, "normal", K, M=M, verbose=False)
    pv = result.stats["pvalue"].values
    assert_allclose(pv[:2], [8.159539103135342e-05, 0.10807353641893498])


def test_qtl_st_scan_lmm_repeat_samples_by_index():
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

    result = st_scan(X, y, "normal", K, M=M, verbose=False)
    pv = result.stats["pvalue"]
    assert_allclose(pv.values[0], 0.9920306566395604)

    ix_best_snp = argmin(array(result.stats["pvalue"]))

    M = concatenate((M, X.loc[:, [ix_best_snp]]), axis=1)
    M = DataFrame(data=M, index=samples)

    result = st_scan(X, y, "normal", K, M=M, verbose=False)
    pv = result.stats["pvalue"]
    assert_allclose(pv[ix_best_snp], 1.0)
    assert_allclose(pv.values[0], 0.6684700834450028)


def test_qtl_st_scan_glmm_binomial():
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

    result = st_scan(X, successes, ("binomial", ntrials), K, verbose=False)
    pv = result.stats["pvalue"]
    assert_allclose(pv, [0.409114, 0.697728], atol=1e-6, rtol=1e-6)


def test_qtl_st_scan_glmm_wrong_dimensions():
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
    st_scan(X, successes, ("binomial", ntrials), K, M=M, verbose=False)


def test_qtl_st_scan_glmm_bernoulli():
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

    result = st_scan(X, successes, "bernoulli", K, verbose=False)
    pv = result.stats["pvalue"]
    assert_allclose(pv, [0.3824950223418756, 0.39213078238911203], atol=1e-5, rtol=1e-5)


def test_qtl_st_scan_glmm_bernoulli_nokinship():
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

    result = st_scan(X, successes, "bernoulli", verbose=False)
    pv = result.stats["pvalue"]
    assert_allclose(pv, [0.9259612341394918, 0.1767987580861164], atol=1e-5, rtol=1e-5)


def test_qtl_st_scan_lm():
    random = RandomState(0)
    nsamples = 50

    G = random.randn(50, 100)

    y = dot(G, random.randn(100)) / sqrt(100) + 0.2 * random.randn(nsamples)

    M = G[:, :5]
    X = G[:, 5:]
    result = st_scan(X, y, "normal", M=M, verbose=False)
    pv = result.stats["pvalue"]
    assert_allclose(pv[:2], [0.133021212899, 0.0312315507648], rtol=1e-4)


def test_qtl_st_scan_gmm_binomial():
    random = RandomState(0)
    nsamples = 50

    X = random.randn(nsamples, 2)
    ntrials = random.randint(1, nsamples, nsamples)
    z = dot(X, random.randn(2))

    successes = zeros(len(ntrials), int)
    for i in range(len(ntrials)):
        for _ in range(ntrials[i]):
            successes[i] += int(z[i] + 0.5 * random.randn() > 0)

    result = st_scan(X, successes, ("binomial", ntrials), verbose=False)
    pv = result.stats["pvalue"]
    assert_allclose(pv, [3.527941e-01, 6.197322e-12], rtol=1e-5, atol=1e-5)


if __name__ == "__main__":
    test_qtl_st_scan_glmm_bernoulli_nokinship()
