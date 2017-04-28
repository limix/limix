import numpy as np
import scipy.stats as st
from numpy import mean as _mean
from numpy import std as _std
from numpy import arange


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


def plot_normal(x, bins=20, nstd=2, figure=None):
    r"""Plot a fit of a normal distribution to the data in x.

    Args:
        x (array_like): values to be fitted.
        bins (int): number of histogram bins.
        nstd (float): standard deviation multiplier.
        figure (:class:`matplotlib.figure.Figure`): matplotlib figure or
                                                    ``None``.

    Returns:
        :class:`matplotlib.figure.Figure`: figure to be shown.

    Examples
    --------

        .. plot::

            from numpy.random import RandomState
            from limix.plot import plot_normal
            from matplotlib import pyplot as plt
            random = RandomState(10)
            x = random.randn(100)
            f = plot_normal(x, nstd=2)
            plt.show(f)
    """
    import matplotlib.pyplot as plt

    if figure is None:
        figure = plt.figure()

    ax = figure.gca()

    mean_x = _mean(x)
    std_x = _std(x)

    xvals = arange(mean_x - 5 * std_x, mean_x + 5 * std_x, .001)
    yvals = st.norm.pdf(xvals, mean_x, std_x)

    ax.hist(x, bins, normed=True)

    ax.plot(xvals, yvals, color='red')

    draw_normal(ax, mean_x, std_x, nstd, 'red')

    return figure
