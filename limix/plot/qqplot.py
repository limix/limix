from __future__ import division

from numpy import sum as npsum
from numpy import (
    arange, ascontiguousarray, flipud, linspace, log10, ones, percentile,
    searchsorted, sort, where
)
from scipy.special import betaincinv


def plot_qqplot(df, alpha=0.05, style=None, ax=None):
    r"""Quantile-Quantile plot of observed p-values versus theoretical ones.

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        Data frame.
    alpha : float
        Significance level defining the band boundary. Defaults to ``0.05``.
    style : dict
        Keyword arguments forwarded to :func:`matplotlib.axes.Axes.plot`
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

        import pandas as pd
        from limix.plot import plot_qqplot
        from numpy.random import RandomState
        from matplotlib import pyplot as plt
        random = RandomState(1)

        pv0 = random.rand(10000)
        pv1 = random.rand(10000)

        data = dict(pv=list(pv0) + list(pv1),
                    label=['label0'] * len(pv0) + ['label1'] * len(pv1))
        plot_qqplot(pd.DataFrame(data=data))
    """

    import matplotlib.pyplot as plt

    ax = plt.gca() if ax is None else ax

    labels = list(df['label'].unique())
    if style is None:
        style = {label: dict() for label in labels}

    for label in labels:
        df0 = df.loc[df['label'] == label, ['pv']]['pv']
        pv = df0.reset_index(drop=True).values
        pv = sort(pv)
        ok = _subsample(pv)

        qnull = -log10((0.5 + arange(len(pv))) / len(pv))
        qemp = -log10(pv)

        ax.plot(
            qnull[ok],
            qemp[ok],
            '-o',
            label=label,
            markeredgewidth=0.5,
            markersize=2.2,
            linewidth=0.65,
            **style.get(label))

    ax.plot([0, qnull.max()], [0, qnull.max()], 'r')
    _plot_confidence_band(ok, qnull, alpha, ax)
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
    ax.legend(loc=2)


def _expected(n):
    lnpv = linspace(1 / (n + 1), n / (n + 1), n, endpoint=True)
    return flipud(-log10(lnpv))


def _rank_confidence_band(nranks, significance_level, ok):
    alpha = significance_level

    k0 = arange(1, nranks + 1)
    k1 = flipud(k0).copy()

    k0 = ascontiguousarray(k0[ok])
    k1 = ascontiguousarray(k1[ok])

    top = betaincinv(k0, k1, 1 - alpha)
    bottom = betaincinv(k0, k1, alpha)

    return (bottom, top)


def _plot_confidence_band(ok, null_qvals, significance_level, ax):

    (bo, to) = _rank_confidence_band(len(null_qvals), significance_level, ok)

    bo = -log10(bo)
    to = -log10(to)

    ax.fill_between(
        null_qvals[ok],
        bo,
        to,
        facecolor='black',
        edgecolor='black',
        alpha=0.15)


def _subsample(pvalues):
    resolution = 500

    if len(pvalues) <= resolution:
        return ones(len(pvalues), dtype=bool)

    ok = pvalues <= percentile(pvalues, 0.1)
    nok = ~ok

    qv = -log10(pvalues[nok])
    qv_min = qv[-1]
    qv_max = qv[0]

    snok = npsum(nok)

    resolution = min(snok, resolution)

    qv_chosen = linspace(qv_min, qv_max, resolution)
    pv_chosen = 10**(-qv_chosen)

    idx = searchsorted(pvalues[nok], pv_chosen)
    ok[where(nok)[0][idx]] = True

    ok[0] = True
    ok[-1] = True

    return ok
