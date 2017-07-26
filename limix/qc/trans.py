from __future__ import division

from limix.util.npy_dask import (
    asarray, isnan, zeros_like, abs, nanmean, nanstd, clip)
from numpy import inf
from numpy_sugar import epsilon
from brent_search import brent


def mean_standardize(X, axis=None, out=None):
    r"""Zero-mean and one-deviation normalisation.

    Normalise in such a way that the mean and variance are equal to zero and
    one.
    This transformation is taken over the flattened array by default, otherwise
    over the specified axis.
    Missing values represented by ``NaN`` are ignored.

    It works well with `Dask`_ array.

    Parameters
    ----------
    X : array_like
        Array to be normalised.
    axis : None or int, optional
        Axis along which the normalisation will take place. The default is to
        normalise the flattened array.

    Returns
    -------
    array_like
        Normalized array.

    Examples
    --------

    .. doctest::

        >>> import limix
        >>> from numpy import arange
        >>>
        >>> X = arange(15).reshape((5, 3))
        >>> print(X)
        [[ 0  1  2]
         [ 3  4  5]
         [ 6  7  8]
         [ 9 10 11]
         [12 13 14]]
        >>> X = limix.qc.mean_standardize(X, axis=0)
        >>> print(X)
        [[-1.4142 -1.4142 -1.4142]
         [-0.7071 -0.7071 -0.7071]
         [ 0.      0.      0.    ]
         [ 0.7071  0.7071  0.7071]
         [ 1.4142  1.4142  1.4142]]
    """

    X = asarray(X, dtype=float)

    if axis is None:
        X = X.ravel()
        axis = 0

    shape = X.shape
    nshape = shape[:axis] + (1,) + shape[axis + 1:]

    X = X - nanmean(X, axis=axis).reshape(nshape)
    d = nanstd(X, axis=axis).reshape(nshape)
    d = clip(d, epsilon.tiny, inf)
    X /= d

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
            \frac{x^{\lambda} - 1}{\lambda}, & \text{if } \lambda > 0; \\
            \log(x), & \text{if } \lambda = 0.
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


if __name__ == '__main__':
    import limix
    import numpy as np
    import dask.array as da

    np.random.seed(0)
    X = np.random.randn(5, 3)

    Y = limix.qc.mean_standardize(X)
    print(Y.mean())
    print(Y.std())

    axis = 0
    Y = limix.qc.mean_standardize(X, axis=axis)
    print(Y.mean(axis))
    print(Y.std(axis))

    axis = 1
    Y = limix.qc.mean_standardize(X, axis=axis)
    print(Y.mean(axis))
    print(Y.std(axis))

    X = da.from_array(X, chunks=(2, 1))
    Y = limix.qc.mean_standardize(X)
    print(Y.mean().compute())
    print(Y.std().compute())

    axis = 0
    Y = limix.qc.mean_standardize(X, axis=axis)
    print(Y.mean(axis).compute())
    print(Y.std(axis).compute())

    axis = 1
    Y = limix.qc.mean_standardize(X, axis=axis)
    print(Y.mean(axis).compute())
    print(Y.std(axis).compute())
