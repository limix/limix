from __future__ import division

from numpy import dot, zeros
from numpy.random import RandomState
from numpy.testing import assert_allclose

from limix.qtl import iscan
from pandas import DataFrame


# pd.set_option('precision', 5)


def test_qtl_interact():
    random = RandomState(0)
    nsamples = 50

    X = random.randn(50, 10)
    G = random.randn(50, 100)
    K = dot(G, G.T)
    ntrials = random.randint(1, 100, nsamples)
    z = dot(G, random.randn(100)) / 10

    successes = zeros(len(ntrials), int)
    for i in range(len(ntrials)):
        for j in range(ntrials[i]):
            successes[i] += int(z[i] + 0.5 * random.randn() > 0)

    y = successes / ntrials

    inter = random.randn(50, 3)

    index = ['sample%02d' % i for i in range(X.shape[0])]
    cols = ['SNP%02d' % i for i in range(X.shape[1])]
    X = DataFrame(data=X, index=index, columns=cols)

    cols = ['inter%02d' % i for i in range(inter.shape[1])]
    inter = DataFrame(data=inter, index=index, columns=cols)

    model = iscan(X, y, 'normal', inter, K, verbose=False)

    print(model.variant_pvalues)

    # pv = lmm.variant_pvalues
    # assert_allclose(pv, [0.44255951309982378, 0.67960798630622032], rtol=1e-3)
    # print(model.variant_pvalues.loc['SNP02', 'inter01'])
    assert_allclose(model.variant_pvalues.loc['SNP02', 'inter01'],
                    [0.72881611358])
