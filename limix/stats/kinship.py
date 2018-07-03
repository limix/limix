from __future__ import division

import numpy as np
from numpy import ascontiguousarray, sqrt, zeros
from tqdm import tqdm


def linear_kinship(G, out=None, verbose=True):
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
        >>> K = linear_kinship(X, verbose=False)
        >>> print(array_str(K, precision=4))
        [[ 0.9131 -0.1928 -0.3413 -0.379 ]
         [-0.1928  0.8989 -0.2356 -0.4704]
         [-0.3413 -0.2356  0.9578 -0.3808]
         [-0.379  -0.4704 -0.3808  1.2302]]
    """
    (n, p) = G.shape
    if out is None:
        out = zeros((n, n))

    chunks = _get_chunks(G)

    start = 0
    for chunk in tqdm(chunks, desc="Kinship", disable=not verbose):
        end = start + chunk

        g = np.asarray(G[:, start:end])
        m = np.nanmean(g, 0)
        g = np.where(np.isnan(g), m, g)
        g = g - m
        g /= np.std(g, 0)
        g /= sqrt(p)

        out += ascontiguousarray(g.dot(g.T), float)

        start = end

    return out


def _get_chunks(G):
    if hasattr(G, "chunks"):
        return G.chunks[1]

    siz = G.shape[1] // 100
    sizl = G.shape[1] - siz * 100
    chunks = [siz] * 100
    if sizl > 0:
        chunks += [sizl]
    return tuple(chunks)
