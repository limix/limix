def quantile_gaussianize(X):
    r"""Normalize a sequence of values via rank and Normal c.d.f.

    Parameters
    ----------
    X : array_like
        Array of values.

    Returns
    -------
    array_like
        Gaussian-normalized values.

    Examples
    --------
    .. doctest::

        >>> from limix.qc import quantile_gaussianize
        >>> from numpy import array_str
        >>>
        >>> qg = quantile_gaussianize([-1, 0, 2])
        >>> print(qg) # doctest: +FLOAT_CMP
        [-0.6744897501960817, 0.0, 0.6744897501960817]
    """
    from .._bits import dask, xarray

    orig_shape = None
    if hasattr(X, "shape"):
        orig_shape = _get_shape(X)
        if X.ndim == 1:
            X = X.reshape(orig_shape + (1,))

    if hasattr(X, "astype"):
        X = X.astype(float)

    if dask.is_array(X):
        X = _quantile_gaussianize_dask(X)
    elif dask.is_data_frame(X):
        x = _quantile_gaussianize_dask(X)
        import dask.dataframe as dd

        x = x.rechunk([(x.shape[0] + 1) // X.index.npartitions, x.chunks[1]])
        X = dd.from_dask_array(x, columns=X.columns, index=X.index)
    elif xarray.is_data_array(X):
        data = X.data

        if dask.is_array(data):
            data = _quantile_gaussianize_dask(data)
        else:
            data = _quantile_gaussianize_ndarray(data)

        X.data = data
    else:
        if hasattr(X, "to_numpy"):
            x = X.to_numpy()
        else:
            x = X

        X[:] = _quantile_gaussianize_ndarray(x)

    if orig_shape is not None and hasattr(X, "reshape"):
        X = X.reshape(orig_shape)

    return X


def _get_shape(x):
    from .._bits import dask

    if dask.is_array(x) or dask.is_data_frame(x):
        import dask.array as da

        return da.compute(*x.shape)
    return x.shape


def _quantile_gaussianize_ndarray(X):
    from scipy.stats import norm
    from numpy import isfinite, asanyarray
    from numpy.ma import masked_invalid
    from bottleneck import nanrankdata
    from numpy import apply_along_axis

    X = asanyarray(X, float)
    orig_shape = X.shape
    if X.ndim == 1:
        X = X.reshape(orig_shape + (1,))
    X = masked_invalid(X)
    X *= -1
    X = nanrankdata(X, axis=0)
    X = X / (isfinite(X).sum(axis=0) + 1)
    X = apply_along_axis(norm.isf, 0, X)
    return X.reshape(orig_shape)


def _quantile_gaussianize_dask(x):
    import dask.array as da
    from scipy.stats import norm
    from bottleneck import nanrankdata
    from .._bits import dask

    if hasattr(x, "to_dask_array"):
        x = x.to_dask_array(lengths=True)

    x = dask.array_shape_reveal(x)
    shape = da.compute(*x.shape)
    x = da.ma.masked_array(x)
    x *= -1
    x = da.apply_along_axis(_dask_apply, 0, x, nanrankdata, shape[0])
    x = x / (da.isfinite(x).sum(axis=0) + 1)
    x = da.apply_along_axis(_dask_apply, 0, x, norm.isf, shape[0])

    return x


def _dask_apply(x, func1d, length):
    from numpy import resize

    try:
        x = func1d(x)
        x = resize(x, length)
    except:
        breakpoint()
    return x
