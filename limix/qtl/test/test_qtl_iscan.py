import pytest
import scipy.stats as st
from limix.qc import normalise_covariance
from limix.qtl import iscan
from limix.stats import linear_kinship
from numpy import (
    argmin,
    array,
    asarray,
    concatenate,
    dot,
    exp,
    eye,
    kron,
    nan,
    ones,
    reshape,
    sqrt,
    zeros,
)
from numpy.random import RandomState
from numpy.testing import assert_allclose, assert_array_equal
from pandas import DataFrame


def _test_qtl_iscan_three_hypothesis(lik):
    random = RandomState(0)
    n = 30
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
    print(r)


def test_qtl_iscan_three_hypothesis():
    _test_qtl_iscan_three_hypothesis("normal")
    _test_qtl_iscan_three_hypothesis("poisson")
    _test_qtl_iscan_three_hypothesis("bernoulli")
    _test_qtl_iscan_three_hypothesis("probit")
    _test_qtl_iscan_three_hypothesis("binomial")


def test_qtl_iscan_two_hypothesis_1vs0():
    random = RandomState(4)
    n = 30
    ncovariates = 3

    M = random.randn(n, ncovariates)

    E0 = random.randint(0, 2, (n, 1)).astype(float)

    G = random.randn(n, 4)

    K = random.randn(n, n + 1)
    K = normalise_covariance(K @ K.T)

    beta = random.randn(ncovariates)
    alpha0 = random.randn(E0.shape[1])

    mvn = random.multivariate_normal
    y = _normalize(M @ beta) + _normalize(E0 @ alpha0)
    y += _normalize(mvn(zeros(n), K + eye(n)))

    idx = [[0, 1], 2, [3]]
    r = iscan(G, y, idx=idx, K=K, M=M, E0=E0, verbose=False)
    print(r)


def test_qtl_iscan_two_hypothesis_2vs0():
    random = RandomState(4)
    n = 30
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
    print(r)


def _normalize(x):
    return (x - x.mean()) / x.std()
