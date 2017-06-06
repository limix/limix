from numpy import ascontiguousarray, copyto, sqrt, zeros
from tqdm import tqdm


def linear_kinship(G, out=None, progress=True):
    r"""Estimate Kinship matrix via linear kernel.

    Examples
    --------
    .. doctest::

        >>> from numpy.random import RandomState
        >>> from limix.stats import linear_kinship
        >>>
        >>> random = RandomState(1)
        >>> X = random.randn(4, 100)
        >>> K = linear_kinship(X, progress=False)
        >>> print(K)
        [[ 27.3944  -5.785  -10.2402 -11.3693]
         [ -5.785   26.9655  -7.068  -14.1125]
         [-10.2402  -7.068   28.7332 -11.425 ]
         [-11.3693 -14.1125 -11.425   36.9068]]
    """
    (n, p) = G.shape
    if out is None:
        out = zeros((n, n))

    nsteps = min(30, p)

    for i in tqdm(range(nsteps), disable=not progress):
        start = i * (p // nsteps)
        stop = min(start + p // nsteps, p)

        G = G - G.mean(0)
        G /= G.std(0)
        G /= sqrt(p)

        out += ascontiguousarray(G.dot(G.T), float)

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
