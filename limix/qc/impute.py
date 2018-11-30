from __future__ import division


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
        X = _impute_dask_array(X)

    elif isinstance(X, xr.DataArray):
        data = X.data

        if isinstance(data, da.Array):
            data = _impute_dask_array(data)
        else:
            data = _impute_ndarray(data)

        X.data = data
    else:
        if hasattr(X, "values"):
            x = X.values
        else:
            x = X
        x = _impute_ndarray(x)

    return X


def _impute_ndarray(x):
    from numpy import isnan, nanmean

    m = nanmean(x, axis=0)
    for i, mi in enumerate(m):
        x[isnan(x[:, i]), i] = mi
    return x


def _impute_dask_array(x):
    import dask.array as da

    m = da.nanmean(x, axis=0).compute()
    start = 0

    arrs = []
    for i in range(len(x.chunks[1])):
        end = start + x.chunks[1][i]
        impute = _get_imputer(m[start:end])
        arrs.append(x[:, start:end].map_blocks(impute, dtype=float))
        start = end
    return da.concatenate(arrs, axis=1)


def _get_imputer(m):
    from numpy import isnan

    def impute(X):
        A = X.copy()

        isn = isnan(A)
        A[:] = 0
        A[isn] = 1

        X[isn] = 0
        X += A * m

        return X

    return impute
