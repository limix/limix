from __future__ import division

from numpy import isnan, nanmean


def _get_imputer(m):
    def impute(X):
        A = X.copy()

        isn = isnan(A)
        A[:] = 0
        A[isn] = 1

        X[isn] = 0
        X += A * m

        return X

    return impute


def mean_impute(G):
    r"""Column-wise impute ``NaN`` values by column mean.

    Parameters
    ----------
    G : array_like
        Bidimensional array to be imputed.

    Returns
    -------
    array_like
        Imputed array.

    Examples
    --------
    .. doctest::

        >>> from numpy.random import RandomState
        >>> from numpy import nan
        >>> from limix.qc import mean_impute
        >>>
        >>> random = RandomState(0)
        >>> X = random.randn(5, 2)
        >>> X[0, 0] = nan
        >>>
        >>> print(mean_impute(X))
        [[ 0.9233  0.4002]
         [ 0.9787  2.2409]
         [ 1.8676 -0.9773]
         [ 0.9501 -0.1514]
         [-0.1032  0.4106]]
    """
    import dask.array as da

    if isinstance(G, da.Array):
        m = da.nanmean(G, axis=0).compute()
        start = 0

        arrs = []
        for i in range(len(G.chunks[1])):
            end = start + G.chunks[1][i]
            impute = _get_imputer(m[start:end])
            arrs.append(G[:, start:end].map_blocks(impute, dtype=float))
            start = end
        G = da.concatenate(arrs, axis=1)
    else:
        m = nanmean(G, axis=0)
        for i in range(len(m)):
            G[isnan(G[:, i]), i] = m[i]

    return G
