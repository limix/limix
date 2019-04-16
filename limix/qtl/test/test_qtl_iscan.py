import scipy.stats as st
from numpy import exp, eye, zeros
from numpy.random import RandomState

from limix.qc import normalise_covariance
from limix.qtl import iscan


def _test_qtl_iscan_three_hypotheses(lik):
    random = RandomState(0)
    n = 10
    ncovariates = 3

    M = random.randn(n, ncovariates)

    E0 = random.randint(0, 2, (n, 1)).astype(float)
    E1 = random.randint(0, 2, (n, 2)).astype(float)

    G = random.randn(n, 4)

    K = random.randn(n, n + 1)
    K = normalise_covariance(K @ K.T)

    beta = random.randn(ncovariates)
    alpha0 = random.randn(E0.shape[1])
    alpha1 = random.randn(E1.shape[1])

    mvn = random.multivariate_normal
    y = _normalize(M @ beta) + _normalize(E0 @ alpha0) + _normalize(E1 @ alpha1)
    y += _normalize(mvn(zeros(n), K + 0.2 * eye(n)))

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

    r = iscan(G, y, lik=lik, idx=idx, K=K, M=M, E0=E0, E1=E1, verbose=False)
    str(r)

    r = iscan(G, y, lik=lik, idx=idx, K=K, M=M, E1=E1, verbose=False)
    str(r)


def test_qtl_iscan_three_hypotheses():
    _test_qtl_iscan_three_hypotheses("normal")
    _test_qtl_iscan_three_hypotheses("poisson")
    _test_qtl_iscan_three_hypotheses("bernoulli")
    _test_qtl_iscan_three_hypotheses("probit")
    _test_qtl_iscan_three_hypotheses("binomial")


def test_qtl_iscan_two_hypotheses():
    random = RandomState(4)
    n = 10
    ncovariates = 3

    M = random.randn(n, ncovariates)

    E1 = random.randint(0, 2, (n, 2)).astype(float)

    G = random.randn(n, 4)

    K = random.randn(n, n + 1)
    K = normalise_covariance(K @ K.T)

    beta = random.randn(ncovariates)
    alpha1 = random.randn(E1.shape[1])

    mvn = random.multivariate_normal
    y = _normalize(M @ beta) + _normalize(E1 @ alpha1)
    y += _normalize(mvn(zeros(n), K + eye(n)))

    idx = [[0, 1], 2, [3]]
    r = iscan(G, y, idx=idx, K=K, M=M, E1=E1, verbose=False)
    str(r)

    r = iscan(G, y, K=K, M=M, E1=E1, verbose=False)
    str(r)


def _normalize(x):
    return (x - x.mean()) / x.std()
