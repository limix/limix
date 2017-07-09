from numpy import ascontiguousarray, copyto, zeros
from tqdm import tqdm


def linear_kinship(G, out=None, verbose=True):
    r"""Estimate Kinship matrix via linear kernel.

    Examples
    --------
    .. doctest::

        >>> from numpy.random import RandomState
        >>> from limix.stats import linear_kinship
        >>>
        >>> random = RandomState(1)
        >>> X = random.randn(4, 100)
        >>> K = linear_kinship(X, verbose=False)
        >>> print(K)
        [[ 0.9131 -0.1928 -0.3413 -0.379 ]
         [-0.1928  0.8989 -0.2356 -0.4704]
         [-0.3413 -0.2356  0.9578 -0.3808]
         [-0.379  -0.4704 -0.3808  1.2302]]
    """
    (n, p) = G.shape
    if out is None:
        out = zeros((n, n))

    nsteps = min(100, p)

    for i in tqdm(range(nsteps), disable=not verbose):
        start = i * (p // nsteps)
        stop = min(start + p // nsteps, p)

        X = G[:, start:stop]
        X = X - X.mean(0)
        X /= X.std(0)

        out += ascontiguousarray(X.dot(X.T), float)

    out /= p

    return out


def gower_norm(K, out=None):
    r"""Perform Gower rescaling of covariance matrix K.

    The rescaled covariance matrix has sample variance of 1.

    Examples
    --------

    .. doctest::

        >>> from numpy.random import RandomState
        >>> from limix.stats import gower_norm
        >>> import scipy as sp
        >>>
        >>> random = RandomState(1)
        >>> X = random.randn(4, 4)
        >>> K = sp.dot(X,X.T)
        >>> Z = random.multivariate_normal(sp.zeros(4), K, 50)
        >>> print("%.3f" % sp.mean(Z.var(1,ddof=1)))
        2.335
        >>> Kn = gower_norm(K)
        >>> Zn = random.multivariate_normal(sp.zeros(4), Kn, 50)
        >>> print("%.3f" % sp.mean(Zn.var(1, ddof=1)))
        0.972
    """

    c = (K.shape[0] - 1) / (K.trace() - K.mean(0).sum())
    if out is None:
        return c * K

    copyto(out, K)
    out *= c
