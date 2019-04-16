from __future__ import division


def linear_kinship(G, out=None, verbose=True):
    """
    Estimate Kinship matrix via linear kernel.

    Let 𝑑 be the number of columns of ``G``. The resulting matrix is given by:

    .. math::

        𝙺 = 𝚇𝚇ᵀ/𝑑

    where

    .. math::

        𝚇ᵢⱼ = (𝙶ᵢⱼ - 𝑚ⱼ) / 𝑠ⱼ

    is the matrix ``G`` column-wise normalized by means 𝑚ⱼ and standard deviations 𝑠ⱼ.
    NaNs are ignored so as to produce matrix ``K`` having only real numbers.

    This functions is specially designed to also handle large matrices ``G`` that would
    otherwise require a large amount of memory if it were to be loaded in memory first.
    For those cases, libraries like Dask come in handy.

    Parameters
    ----------
    G : array_like
        Samples-by-variants matrix.
    out : ndarray
        A location into which the result is stored.
    verbose : bool, optional
        ``True`` for showing progress; ``False`` otherwise. Defauts to ``True``.

    Examples
    --------
    .. doctest::

        >>> from numpy.random import RandomState
        >>> from limix.stats import linear_kinship
        >>>
        >>> random = RandomState(1)
        >>> X = random.randn(4, 100)
        >>> K = linear_kinship(X, verbose=False)
        >>> print(K) # doctest: +FLOAT_CMP
        [[ 0.91314823 -0.19283362 -0.34133897 -0.37897564]
         [-0.19283362  0.89885153 -0.2356003  -0.47041761]
         [-0.34133897 -0.2356003   0.95777313 -0.38083386]
         [-0.37897564 -0.47041761 -0.38083386  1.23022711]]
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
