import numpy as np
import scipy.stats as st
from numpy import mean as _mean
from numpy import std as _std
from numpy import arange


def plot_normal(x, bins=20, nstd=2, style=None, ax=None):
    r"""Plot a fit of a normal distribution to the data in x.

    Parameters
    ----------
    x : array_like
        Values to be fitted.
    bins : int
        Number of histogram bins.
    nstd : float)
        Standard deviation multiplier.
    style : dict
        Keyword arguments forwarded to the :func:`matplotlib.axes.Axes.plot`
        function.
    ax : :class:`matplotlib.axes.Axes`
        The target handle for this figure. If ``None``, the current axes is
        set.

    Returns
    -------
    :class:`matplotlib.axes.Axes`
        Axes.

    Examples
    --------
    .. plot::

        from numpy.random import RandomState
        from matplotlib import pyplot as plt
        from limix.plot import plot_normal

        random = RandomState(10)
        x = random.randn(100)
        ax = plot_normal(x, nstd=2)
        plt.tight_layout()
        plt.show()
    """
    import matplotlib.pyplot as plt

    ax = plt.gca() if ax is None else ax

    if style is None:
        style = dict()

    if 'color' not in style:
        style['color'] = 'red'

    mean_x = _mean(x)
    std_x = _std(x)

    xvals = arange(mean_x - 5 * std_x, mean_x + 5 * std_x, .001)
    yvals = st.norm.pdf(xvals, mean_x, std_x)

    ax.hist(x, bins, normed=True)

    ax.plot(xvals, yvals, **style)

    draw_normal(ax, mean_x, std_x, nstd, 'red')

    return ax


def draw_normal(axis, mean, scale, nstd, color):
    max_pdf = st.norm.pdf(mean, mean, scale)

    axis.plot([mean, mean], [0, max_pdf], color=color, linestyle="--")

    axis.annotate(
        '$\mu$',
        xy=(mean + 0.6 * scale, max_pdf),
        horizontalalignment='center',
        verticalalignment='bottom',
        fontsize=15,
        color=color)

    top = st.norm.pdf(mean + nstd * scale, mean, scale)
    left = mean - nstd * scale
    right = mean + nstd * scale

    axis.plot([right, right], [0, top], color=color, linestyle="--")

    axis.plot([left, left], [0, top], color=color, linestyle="--")

    if int(nstd) == nstd:
        mu_sigma = '$\mu+%d\sigma$' % nstd
    else:
        mu_sigma = '$\mu+%.1f\sigma$' % nstd

    axis.annotate(
        mu_sigma,
        xy=(mean + (1.2 + nstd) * scale, top),
        horizontalalignment='center',
        verticalalignment='bottom',
        fontsize=15,
        color=color)
