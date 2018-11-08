from __future__ import division

from numpy import isnan, nanmean


def mean_impute(X):
    r"""Column-wise impute ``NaN`` values by column mean.

    It works well with `Dask`_ array.

    Parameters
    ----------
    X : array_like
        Matrix to have its missing values imputed.

    Returns
    -------
    array_like
        Imputed array.

    Examples
    --------
    .. doctest::

        >>> from numpy.random import RandomState
        >>> from numpy import nan, array_str
        >>> from limix.qc import mean_impute
        >>>
        >>> random = RandomState(0)
        >>> X = random.randn(5, 2)
        >>> X[0, 0] = nan
        >>>
        >>> print(array_str(mean_impute(X), precision=4))
        [[ 0.9233  0.4002]
         [ 0.9787  2.2409]
         [ 1.8676 -0.9773]
         [ 0.9501 -0.1514]
         [-0.1032  0.4106]]

    .. _Dask: https://dask.pydata.org/
    """
    import dask.array as da
    import xarray as xr

    if isinstance(X, da.Array):
        m = da.nanmean(X, axis=0).compute()
        start = 0

        arrs = []
        for i in range(len(X.chunks[1])):
            end = start + X.chunks[1][i]
            impute = _get_imputer(m[start:end])
            arrs.append(X[:, start:end].map_blocks(impute, dtype=float))
            start = end
        X = da.concatenate(arrs, axis=1)
    elif isinstance(X, xr.DataArray):
        data = X.data
        m = da.nanmean(data, axis=0).compute()
        start = 0

        arrs = []
        for i in range(len(data.chunks[1])):
            end = start + data.chunks[1][i]
            impute = _get_imputer(m[start:end])
            arrs.append(data[:, start:end].map_blocks(impute, dtype=float))
            start = end
        data = da.concatenate(arrs, axis=1)
        X.data = data
    else:
        if hasattr(X, "values"):
            x = X.values
        else:
            x = X
        m = nanmean(x, axis=0)
        for i, mi in enumerate(m):
            x[isnan(x[:, i]), i] = mi

    return X


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
