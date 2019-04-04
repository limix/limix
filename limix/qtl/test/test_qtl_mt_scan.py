from limix.qtl import mt_scan
from numpy import reshape, kron, eye
from numpy import concatenate
from numpy.random import RandomState
from numpy.testing import assert_allclose
import scipy.stats as st
from limix.qc import normalise_covariance


def test_qtl_mt_scan():
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
    r = mt_scan(G, Y, idx=idx, K=K, M=M, A=A, A0=A0, A1=A1, verbose=False)
    # breakpoint()
    # print()
    # print(r.stats)
    # print(r.effsizes["h1"])
    # print(r.effsizes["h2"])
    # print(r.h0.effsizes)
    # print(r.h0.variances)
    # print(r.h0.lml)
    # print(r.null_beta())

    # A = random.randn(3, 3)
    # A = A @ A.T
    # M = random.randn(n, 2)


def vec(x):
    return reshape(x, (-1,) + x.shape[2:], order="F")


def unvec(x, shape):
    return reshape(x, shape, order="F")
