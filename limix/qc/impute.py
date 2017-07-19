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
        >>> from limix.qc import boxcox
        >>>
        >>> random = RandomState(0)
        >>> X = random.randn(5, 2)
        >>> X[0, 0] = nan
        >>>
        >>> print(boxcox(X))
        [[    nan  0.9544]
         [ 1.0986  1.6946]
         [ 2.0183  0.    ]
         [ 1.0691  0.644 ]
         [ 0.      0.9597]]
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
        m = nanmean(G, axis=0).compute()
        for i in range(len(m)):
            G[isnan(G[:, i]), i] = m[i]

    return G
