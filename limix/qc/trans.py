from __future__ import division

from numpy import asarray, isnan, zeros_like


def mean_standardize(X, inplace=False):
    r"""Zero-mean and one-deviation normalization.

    Standardize Y in a way that is robust to missing values.

    Parameters
    ----------
    Y : array_like
        Array to be normalized.
    inplace : bool
        Whether to operate in-place.

    Returns
    -------
    array_like
        Normalized array.
    """
    import dask.array as da

    if isinstance(X, da.Array):
        X = X.astype(float)
        X = X - da.nanmean(X, axis=0)
        X = X / da.nanstd(X, axis=0)
    else:
        X = asarray(X, float)
        if inplace:
            X -= X.nanmean(0)
        else:
            X = X - X.nanmean(0)
        dev = X.nanstd(0)
        X[:, dev > 0] /= dev[dev > 0]

    return X


def quantile_gaussianize(x):
    r"""Normalize a sequence of values via rank and Normal c.d.f.

    Parameters
    ----------
    x : array_like
        Sequence of values.

    Returns
    -------
    array_like
        Gaussian-normalized values.

    Examples
    --------
    .. doctest::

        >>> from limix.qc import quantile_gaussianize
        >>> print(quantile_gaussianize([-1, 0, 2]))
        [-0.6745  0.      0.6745]
    """

    from scipy_sugar.stats import quantile_gaussianize as qg

    return qg(x)


def boxcox(X):
    r"""Gaussianize X using the Box-Cox transformation.

    Each phentoype is brought to a positive schale by first subtracting the
    minimum value and adding 1.
    Then each phenotype is transformed by the Box-Cox transformation.

    Parameters
    ----------
    X : array_like
        Samples by phenotypes.

    Returns
    -------
    array_like
        Box-Cox power transformed array.

    Examples
    --------
    .. doctest::

        >>> from numpy.random import RandomState
        >>> from limix.qc import boxcox
        >>>
        >>> random = RandomState(0)
        >>> X = random.randn(5, 2)
        >>>
        >>> print(boxcox(X))
        [[ 2.7136  0.9544]
         [ 1.3844  1.6946]
         [ 2.9066  0.    ]
         [ 1.3407  0.644 ]
         [ 0.      0.9597]]
    """
    from scipy.stats import boxcox as sp_boxcox

    X_transformed = zeros_like(X)
    for i in range(X.shape[1]):
        i_nan = isnan(X[:, i])
        values = X[~i_nan, i]
        X_transformed[i_nan, i] = X[i_nan, i]
        X_transformed[~i_nan, i], _ = sp_boxcox(values - values.min() + 1.0)
    return X_transformed
