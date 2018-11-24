from __future__ import division


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
    from numpy import sqrt, zeros, asfortranarray, where, asarray, nanmean, std, isnan
    from scipy.linalg.blas import get_blas_funcs
    from tqdm import tqdm

    (n, p) = G.shape
    if out is None:
        out = zeros((n, n), order="F")
    else:
        out = asfortranarray(out)

    chunks = _get_chunks(G)
    gemm = get_blas_funcs("gemm", [out])

    start = 0
    for chunk in tqdm(chunks, desc="Kinship", disable=not verbose):
        end = start + chunk

        g = asarray(G[:, start:end])
        m = nanmean(g, 0)
        g = where(isnan(g), m, g)
        g = g - m
        g /= std(g, 0)
        g /= sqrt(p)

        gemm(1.0, g, g, 1.0, out, 0, 1, 1)

        start = end

    return out


def _get_chunks(G):
    from numpy import isfinite

    if hasattr(G, "chunks") and G.chunks is not None:
        if len(G.chunks) > 1 and all(isfinite(G.chunks[0])):
            return G.chunks[1]

    siz = G.shape[1] // 100
    sizl = G.shape[1] - siz * 100
    chunks = [siz] * 100
    if sizl > 0:
        chunks += [sizl]
    return tuple(chunks)
