from numpy import concatenate, dot, sqrt
from numpy.random import RandomState
from numpy.testing import assert_allclose
from pandas import DataFrame

from limix.qtl import st_scan
from limix.stats import linear_kinship


def test_qtl_lmm():
    random = RandomState(0)
    nsamples = 50

    G = random.randn(50, 100)
    K = linear_kinship(G[:, 0:80], verbose=False)

    y = dot(G, random.randn(100)) / sqrt(100) + 0.2 * random.randn(nsamples)

    M = G[:, :5]
    X = G[:, 68:70]

    model = st_scan(X, y, "normal", K, M=M, verbose=False)
    pv = model.variant_pvalues

    ix_best_snp = pv.argmin().item()

    M = concatenate((M, X[:, [ix_best_snp]]), axis=1)

    model = st_scan(X, y, "normal", K, M=M, verbose=False)
    pv = model.variant_pvalues
    assert_allclose(pv[ix_best_snp], 1.0)


def test_qtl_lmm_nokinship():
    random = RandomState(0)
    nsamples = 50

    G = random.randn(50, 100)
    K = linear_kinship(G[:, 0:80], verbose=False)

    y = dot(G, random.randn(100)) / sqrt(100) + 0.2 * random.randn(nsamples)

    M = G[:, :5]
    X = G[:, 68:70]

    model = st_scan(X, y, "normal", K, M=M, verbose=False)
    pv = model.variant_pvalues.values
    assert_allclose(pv[:2], [8.159539103135342e-05, 0.10807353641893498])


def test_qtl_lmm_repeat_samples_by_index():
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

    model = st_scan(X, y, "normal", K, M=M, verbose=False)
    pv = model.variant_pvalues
    assert_allclose(pv.values[0], 0.9920306566395604)

    ix_best_snp = pv.argmin().item()

    M = concatenate((M, X.loc[:, [ix_best_snp]]), axis=1)
    M = DataFrame(data=M, index=samples)

    model = st_scan(X, y, "normal", K, M=M, verbose=False)
    pv = model.variant_pvalues
    assert_allclose(pv[ix_best_snp], 1.0)
    assert_allclose(pv.values[0], 0.6684700834450028)
