def normalise_covariance(K, out=None):
    """
    Variance rescaling of covariance matrix 𝙺.

    Let n be the number of rows (or columns) of 𝙺 and let
    mᵢ be the average of the values in the i-th column.
    Gower rescaling is defined as

    .. math::

        𝙺(n - 1)/(𝚝𝚛𝚊𝚌𝚎(𝙺) - ∑mᵢ).

    Notes
    -----
    The reasoning of the scaling is as follows.
    Let 𝐠 be a vector of n independent samples and let 𝙲 be the Gower's centering
    matrix.
    The unbiased variance estimator is

    .. math::

        v = ∑ (gᵢ-ḡ)²/(n-1) = 𝚝𝚛𝚊𝚌𝚎((𝐠-ḡ𝟏)ᵀ(𝐠-ḡ𝟏))/(n-1) = 𝚝𝚛𝚊𝚌𝚎(𝙲𝐠𝐠ᵀ𝙲)/(n-1)

    Let 𝙺 be the covariance matrix of 𝐠.
    The expectation of the unbiased variance estimator is

    .. math::

        𝐄[v] = 𝚝𝚛𝚊𝚌𝚎(𝙲𝐄[𝐠𝐠ᵀ]𝙲)/(n-1) = 𝚝𝚛𝚊𝚌𝚎(𝙲𝙺𝙲)/(n-1),

    assuming that 𝐄[gᵢ]=0.
    We thus divide 𝙺 by 𝐄[v] to achieve an unbiased normalisation on the random variable
    gᵢ.

    Parameters
    ----------
    K : array_like
        Covariance matrix to be normalised.
    out : array_like, optional
        Result destination. Defaults to ``None``.

    Examples
    --------
    .. doctest::

        >>> from numpy import dot, mean, zeros
        >>> from numpy.random import RandomState
        >>> from limix.qc import normalise_covariance
        >>>
        >>> random = RandomState(0)
        >>> X = random.randn(10, 10)
        >>> K = dot(X, X.T)
        >>> Z = random.multivariate_normal(zeros(10), K, 500)
        >>> print("%.3f" % mean(Z.var(1, ddof=1)))
        9.824
        >>> Kn = normalise_covariance(K)
        >>> Zn = random.multivariate_normal(zeros(10), Kn, 500)
        >>> print("%.3f" % mean(Zn.var(1, ddof=1)))
        1.008

    .. _Dask: https://dask.pydata.org/
    """
    from numpy import asarray
    import dask.array as da
    from pandas import DataFrame
    import xarray as xr

    if isinstance(K, DataFrame):
        K = K.astype(float)
        trace = K.values.trace()
    elif isinstance(K, da.Array):
        trace = da.diag(K).sum()
    elif isinstance(K, xr.DataArray):
        trace = da.diag(K.data).sum()
    else:
        K = asarray(K, float)
        trace = K.trace()

    c = asarray((K.shape[0] - 1) / (trace - K.mean(axis=0).sum()), float)
    if out is None:
        return K * c

    _copyto(out, K)
    _inplace_mult(out, c)

    return out


def _copyto(dst, src):
    from numpy import copyto, ndarray
    import dask.array as da
    from pandas import DataFrame

    if isinstance(dst, DataFrame):
        copyto(dst.values, src)
    elif isinstance(dst, ndarray) and isinstance(src, da.Array):
        copyto(dst, src.compute())
    else:
        copyto(dst, src)


def _inplace_mult(out, c):
    import dask.array as da

    if isinstance(c, da.Array):
        c = c.compute()
    out *= c
