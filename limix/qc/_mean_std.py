import warnings

from limix._bits import dask, numpy, pandas, xarray


def mean_standardize(X, axis=-1, inplace=False):
    """
    Zero-mean and one-deviation normalisation.

    Normalise in such a way that the mean and variance are equal to zero and
    one. This transformation is taken over the flattened array by default, otherwise
    over the specified axis. Missing values represented by ``NaN`` are ignored.

    Parameters
    ----------
    X : array_like
        Array of values.
    axis : int, optional
        Axis value. Defaults to ``1``.
    inplace : bool, optional
        Defaults to ``False``.

    Returns
    -------
    X : ndarray
        Normalized array.

    Example
    -------
    .. doctest::

        >>> import limix
        >>> from numpy import arange
        >>>
        >>> X = arange(15).reshape((5, 3)).astype(float)
        >>> print(X)
        [[ 0.  1.  2.]
         [ 3.  4.  5.]
         [ 6.  7.  8.]
         [ 9. 10. 11.]
         [12. 13. 14.]]
        >>> X = arange(6).reshape((2, 3)).astype(float)
        >>> X = limix.qc.mean_standardize(X, axis=0)
        >>> print(X) # doctest: +FLOAT_CMP
        [[-1.22474487  0.          1.22474487]
         [-1.22474487  0.          1.22474487]]
    """

    from numpy import issubdtype, integer, asarray

    if hasattr(X, "dtype") and issubdtype(X.dtype, integer):
        raise ValueError("Integer type is not supported.")

    if isinstance(X, (tuple, list)):
        if inplace:
            raise ValueError("Can't use `inplace=True` for {}.".format(type(X)))
        X = asarray(X, float)

    if numpy.is_array(X):
        X = _mean_standardize_numpy(X, axis, inplace)
    elif pandas.is_series(X):
        X = _mean_standardize_pandas_series(X, axis, inplace)
    elif pandas.is_dataframe(X):
        X = _mean_standardize_pandas_dataframe(X, axis, inplace)
    elif dask.is_array(X):
        X = _mean_standardize_dask_array(X, axis, inplace)
    elif dask.is_series(X):
        raise NotImplementedError()
    elif dask.is_dataframe(X):
        X = _mean_standardize_dask_dataframe(X, axis, inplace)
    elif xarray.is_dataarray(X):
        X = _mean_standardize_xarray_dataarray(X, axis, inplace)
    else:
        raise NotImplementedError()

    return X


def _mean_standardize_numpy(X, axis, inplace):
    from numpy import nanmean, nanstd
    from numpy_sugar import epsilon
    from numpy import inf, clip

    orig_shape = X.shape
    if X.ndim == 1:
        X = X.reshape(orig_shape + (1,))

    if not inplace:
        X = X.copy()

    X = X.swapaxes(1, axis)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        X -= nanmean(X, axis=0)
        X /= clip(nanstd(X, axis=0), epsilon.tiny, inf)

    X = X.swapaxes(1, axis)

    return X.reshape(orig_shape)


def _mean_standardize_pandas_series(X, axis, inplace):
    if not inplace:
        X = X.copy()

    a = X.to_numpy()
    _mean_standardize_numpy(a, axis, True)
    X[:] = a

    return X


def _mean_standardize_pandas_dataframe(x, axis, inplace):
    if not inplace:
        x = x.copy()

    a = x.to_numpy()
    _mean_standardize_numpy(a, axis, True)
    x[:] = a

    return x


def _mean_standardize_dask_array(x, axis, inplace):
    import dask.array as da
    from numpy_sugar import epsilon
    from numpy import nanmean, clip, nanstd, inf

    if inplace:
        raise NotImplementedError()

    x = x.swapaxes(1, axis)

    x = dask.array_shape_reveal(x)
    shape = da.compute(*x.shape)

    def func(a):
        a -= nanmean(a, axis=0)
        a /= clip(nanstd(a, axis=0), epsilon.tiny, inf)
        return a

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        x = da.apply_along_axis(_dask_apply, 0, x, func, shape[0])

    return x.swapaxes(1, axis)


def _mean_standardize_dask_dataframe(x, axis, inplace):
    if inplace:
        raise NotImplementedError()

    d = x.to_dask_array(lengths=True)
    orig_chunks = d.chunks
    d = _mean_standardize_dask_array(d, axis, False).rechunk(orig_chunks)
    return d.to_dask_dataframe(columns=x.columns, index=x.index)


def _mean_standardize_xarray_dataarray(X, axis, inplace):
    if not inplace:
        X = X.copy(deep=True)

    data = X.data

    if dask.is_array(data):
        data = _mean_standardize_dask_array(data, axis, inplace)
    else:
        data = _mean_standardize_numpy(data, axis, inplace)

    X.data = data
    return X


def _dask_apply(x, func1d, length):
    from numpy import resize

    x = func1d(x)
    return resize(x, length)
