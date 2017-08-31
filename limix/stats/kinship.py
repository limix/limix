from __future__ import division

from numpy import ascontiguousarray, sqrt, zeros
from tqdm import tqdm


def linear_kinship(G, out=None, progress=True):
    r"""Estimate Kinship matrix via linear kernel.

    Examples
    --------
    .. doctest::

        >>> from numpy.random import RandomState
        >>> from numpy import array_str
        >>> from limix.stats import linear_kinship
        >>>
        >>> random = RandomState(1)
        >>> X = random.randn(4, 100)
        >>> K = linear_kinship(X, progress=False)
        >>> print(array_str(K, precision=4))
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
