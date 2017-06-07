from __future__ import division

from numpy import append, arange, asarray, flipud, linspace, log10, sort
from scipy.special import betaincinv


def qqplot(pv, label='unknown', alphaLevel=0.05, style=dict(), ax=None):
    r"""Produces a Quantile-Quantile plot of the observed P value
        distribution against the theoretical one under the null.

    Parameters
    ----------

    pv : array_like
        P-values.
    distr : {'log10', 'chi2'}
        Scale of the distribution. If 'log10' is specified, the distribution
        of the -log10 P values is considered.
        If the distribution of the corresponding chi2-distributed test
        statistics is considered. Defaults to 'log10'.
    alphaLevel : float
        Significance bound.
    ax : :class:`matplotlib.axes.AxesSubplot`
        The target handle for this figure. If None, the current axes is set.

    Returns
    -------
    :class:`matplotlib.axes.AxesSubplot`
        Axes.

    Examples
    --------

    .. plot::

        from limix.plot import qqplot
        from numpy.random import RandomState
        from matplotlib import pyplot as plt
        random = RandomState(1)

        pv = random.rand(10000)

        fig = plt.figure(1, figsize=(5,5))
        plt.subplot(111)
        qqplot(pv)
        plt.tight_layout()
        plt.show()
    """

    import matplotlib.pyplot as plt

    ax = plt.gca() if ax is None else ax

    pv = asarray(pv, float)

    qnull = -log10((0.5 + arange(len(pv))) / len(pv))
    qemp = -log10(sort(pv))

    ax.plot(qnull, qemp, '.', label=label, **style)
    ax.plot([0, qnull.max()], [0, qnull.max()], 'r')

    _plot_confidence_band(qnull, alphaLevel, ax)
    _set_axis_labels(ax)

    _set_frame(ax)

    return ax


def _set_frame(ax):
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')


def _set_axis_labels(ax):
    ax.set_ylabel('-log$_{10}$pv observed')
    ax.set_xlabel('-log$_{10}$pv expected')


def _expected(n):
    lnpv = linspace(1 / (n + 1), n / (n + 1), n, endpoint=True)
    return flipud(-log10(lnpv))


def _rank_confidence_band(nranks, significance_level):
    alpha = significance_level

    k0 = arange(1, nranks + 1)
    k1 = flipud(k0).copy()

    top = betaincinv(k0, k1, 1 - alpha)
    mean = k0 / (nranks + 1)
    bottom = betaincinv(k0, k1, alpha)

    return (bottom, mean, top)


def _plot_confidence_band(null_pvals, significance_level, ax):

    (bo, me, to) = _rank_confidence_band(len(null_pvals), significance_level)

    me = -log10(me)
    bo = -log10(bo)
    to = -log10(to)

    ax.fill_between(
        null_pvals, bo, to, facecolor='black', edgecolor='black', alpha=0.15)
