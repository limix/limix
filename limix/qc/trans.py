from __future__ import division

from numpy import asarray, isnan, zeros_like, abs
from numpy_sugar import epsilon
from brent_search import brent


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


def boxcox(x):
    r"""Box Cox transformation for normality conformance.

    It applies the power transformation

    .. math::

        f(x) = \begin{cases}
            \frac{x^{\lambda} - 1}{\lambda}, & \text{if} \lambda > 0; \\
            \log(x), & \text{if} \lambda = 0.
        \end{cases}

    to the provided data, hopefully making it more normal distribution-like.
    The :math:`\lambda` parameter is fit by maximum likelihood estimation.

    Parameters
    ----------
    X : array_like
        Data to be transformed.

    Returns
    -------
    array_like
        Box Cox transformed data.

    Examples
    --------
    .. plot::

        import limix
        import numpy as np
        import scipy.stats as stats

        np.random.seed(0)

        x = stats.loggamma.rvs(0.1, size=100)
        y = limix.qc.boxcox(x)

        fig = limix.plot.figure()

        ax1 = fig.add_subplot(211)
        limix.plot.plot_qqnormal(x, ax=ax1)

        ax2 = fig.add_subplot(212)
        limix.plot.plot_qqnormal(y, ax=ax2)

        limix.plot.show()
    """
    from scipy.stats import boxcox_llf
    from scipy.special import boxcox as bc

    x = asarray(x, float)

    m = x.min()
    if m <= 0:
        m = max(abs(m), epsilon.small)
        x = x + m + m / 2

    lmb = brent(lambda lmb: -boxcox_llf(lmb, x), -5, +5)[0]
    return bc(x, lmb)
