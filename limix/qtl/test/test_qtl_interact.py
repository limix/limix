from __future__ import division

from numpy import dot, zeros
from numpy.random import RandomState
from numpy.testing import assert_allclose

from limix.qtl import iscan
from pandas import DataFrame


def test_qtl_interact():
    random = RandomState(0)
    nsamples = 50

    X = random.randn(nsamples, 10)
    G = random.randn(nsamples, 100)
    K = dot(G, G.T)
    ntrials = random.randint(1, 100, nsamples)
    z = dot(G, random.randn(100)) / 10

    successes = zeros(len(ntrials), int)
    for i in range(len(ntrials)):
        for j in range(ntrials[i]):
            successes[i] += int(z[i] + 0.5 * random.randn() > 0)

    y = successes / ntrials

    inter = random.randn(nsamples, 3)

    index = ['sample%02d' % i for i in range(X.shape[0])]
    cols = ['SNP%02d' % i for i in range(X.shape[1])]
    X = DataFrame(data=X, index=index, columns=cols)

    cols = ['inter%02d' % i for i in range(inter.shape[1])]
    inter = DataFrame(data=inter, index=index, columns=cols)

    model = iscan(X, y, 'normal', inter, K, verbose=False)

    print(model.variant_pvalues)

    assert_allclose(model.variant_pvalues.loc['SNP02', 'inter01'],
                    [0.72881611358])
