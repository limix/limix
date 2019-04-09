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


def test_qtl_iscan_three_hypothesis():
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

    mvn = st.multivariate_normal
    y = _normalize(M @ beta) + _normalize(E0 @ alpha0) + _normalize(E1 @ alpha1)
    y += _normalize(mvn(zeros(n), K + eye(n)).rvs())

    idx = [[0, 1], 2, [3]]
    r = iscan(G, y, idx=idx, K=K, M=M, E0=E0, E1=E1, verbose=False)
    print()
    print(r)


def _normalize(x):
    return (x - x.mean()) / x.std()
