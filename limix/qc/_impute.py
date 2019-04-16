import warnings

from limix._bits import dask, get_shape, numpy, pandas, xarray


def mean_impute(X, axis=-1, inplace=False):
    """
    Impute ``NaN`` values.

    It defaults to column-wise imputation.

    Parameters
    ----------
    X : array_like
        Matrix to have its missing values imputed.
    axis : int, optional
        Axis value. Defaults to `-1`.
    inplace : bool, optional
        Defaults to `False`.

    Returns
    -------
    ndarray
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
    from numpy import asarray

    if isinstance(X, (tuple, list)):
        if inplace:
            raise ValueError("Can't use `inplace=True` for {}.".format(type(X)))
        X = asarray(X, float)

    if numpy.is_array(X):
        X = _impute_numpy(X, axis, inplace)
    elif pandas.is_series(X):
        X = _impute_pandas_series(X, axis, inplace)
    elif pandas.is_dataframe(X):
        X = _impute_pandas_dataframe(X, axis, inplace)
    elif dask.is_array(X):
        X = _impute_dask_array(X, axis, inplace)
    elif dask.is_series(X):
        raise NotImplementedError()
    elif dask.is_dataframe(X):
        X = _impute_dask_dataframe(X, axis, inplace)
    elif xarray.is_dataarray(X):
        X = _impute_xarray_dataarray(X, axis, inplace)
    else:
        raise NotImplementedError()

    return X


def _impute_numpy(X, axis, inplace):
    from numpy import isnan, nanmean

    orig_shape = get_shape(X)
    if X.ndim == 1:
        X = X.reshape(orig_shape + (1,))

    if not inplace:
        X = X.copy()

    X = X.swapaxes(-1, axis)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        m = nanmean(X, axis=0)
    for i, mi in enumerate(m):
        X[isnan(X[:, i]), i] = mi

    X = X.swapaxes(-1, axis)

    return X.reshape(orig_shape)


def _impute_pandas_series(x, axis, inplace):
    if not inplace:
        x = x.copy()

    a = x.to_numpy()
    _impute_numpy(a, axis, True)
    x[:] = a

    return x


def _impute_dask_dataframe(x, axis, inplace):
    import dask.dataframe as dd

    if inplace:
        raise NotImplementedError()

    d = x.to_dask_array(lengths=True)
    d = _impute_dask_array(d, axis, False)
    x = dd.from_dask_array(d, columns=x.columns, index=x.index)
    return x


def _impute_pandas_dataframe(x, axis, inplace):
    if not inplace:
        x = x.copy()

    a = x.to_numpy()
    _impute_numpy(a, axis, True)
    x[:] = a

    return x


def _impute_dask_array(x, axis, inplace):
    import dask.array as da

    x = x.swapaxes(-1, axis)
    m = da.nanmean(x, axis=0).compute()
    start = 0

    arrs = []
    for i in range(len(x.chunks[1])):
        end = start + x.chunks[1][i]
        impute = _get_imputer(m[start:end], inplace)
        arrs.append(x[:, start:end].map_blocks(impute, dtype=float))
        start = end
    return da.concatenate(arrs, axis=1).swapaxes(-1, axis)


def _impute_xarray_dataarray(X, axis, inplace):
    if not inplace:
        X = X.copy(deep=True)

    data = X.data

    if dask.is_array(data):
        data = _impute_dask_array(data, axis, inplace)
    else:
        data = _impute_numpy(data, axis, inplace)

    X.data = data
    return X


def _get_imputer(m, inplace):
    from numpy import isnan

    def impute(X):
        A = X.copy()

        isn = isnan(A)
        A[:] = 0
        A[isn] = 1

        if not inplace:
            X = X.copy()

        X[isn] = 0
        X += A * m

        return X

    return impute
