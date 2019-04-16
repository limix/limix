from __future__ import division

from brent_search import brent


def boxcox(x):
    r"""Box-Cox transformation for normality conformance.

    It applies the power transformation

    .. math::

        f(x) = \begin{cases}
            \frac{x^{\lambda} - 1}{\lambda}, & \text{if } \lambda > 0; \\
            \log(x), & \text{if } \lambda = 0.
        \end{cases}

    to the provided data, hopefully making it more normal distribution-like.
    The λ parameter is fit by maximum likelihood estimation.

    Parameters
    ----------
    X : array_like
        Data to be transformed.

    Returns
    -------
    boxcox : ndarray
        Box-Cox transformed data.

    Examples
    --------
    .. plot::

        >>> import limix
        >>> import numpy as np
        >>> import scipy.stats as stats
        ...
        >>> np.random.seed(0)
        ...
        >>> x = stats.loggamma.rvs(0.1, size=100)
        >>> y = limix.qc.boxcox(x)
        ...
        >>> plt = limix.plot.get_pyplot()
        ...
        >>> _, (ax1, ax2) = plt.subplots(2, 1)
        >>> _ = stats.probplot(x, dist=stats.norm, plot=ax1)
        >>> _ = stats.probplot(y, dist=stats.norm, plot=ax2)
        >>> plt.tight_layout()
    """
    import dask.array as da
    import numpy as np

    if isinstance(x, da.Array):
        return _boxcox(da, x)
    return _boxcox(np, x)


def _boxcox(lib, x):
    from numpy_sugar import epsilon
    from scipy.stats import boxcox_llf
    from scipy.special import boxcox as bc

    x = lib.asarray(x).astype(float)

    m = x.min()
    if m <= 0:
        m = max(lib.abs(m), epsilon.small)
        x = x + m + m / 2

    lmb = brent(lambda lmb: -boxcox_llf(lmb, x), -5, +5)[0]
    return bc(x, lmb)
