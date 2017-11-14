from numpy import dot, sqrt, concatenate
from numpy.random import RandomState
from numpy.testing import assert_allclose

from limix.qtl import scan
from limix.stats import linear_kinship


def test_qtl_lmm():
    random = RandomState(0)
    nsamples = 50

    G = random.randn(50, 100)
    K = linear_kinship(G[:, 0:80], verbose=False)

    y = dot(G, random.randn(100)) / sqrt(100) + 0.2 * random.randn(nsamples)

    M = G[:, :5]
    X = G[:, 68:70]

    model = scan(X, y, 'normal', K, M=M, verbose=False)
    pv = model.variant_pvalues
    ix_best_snp = pv.idxmin()
    pos_best_snp = pv.values.argmin()

    M = concatenate((M, X[:, [pos_best_snp]]), axis=1)

    model = scan(X, y, 'normal', K, M=M, verbose=False)
    pv = model.variant_pvalues
    assert_allclose(pv.loc[ix_best_snp], 1.0)
