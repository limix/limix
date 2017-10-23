from __future__ import division

from numpy import sum as npsum
from numpy import (arange, ascontiguousarray, flipud, linspace, log10, ones,
                   percentile, searchsorted, sort, where, atleast_2d)
from scipy.special import betaincinv


def qqplot(df, alpha, style=None, ax=None):
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
    """
    df = _normalise_data(df)

    if style is None:
        style = dict()

    for label, df0 in df.groupby('label'):
        pv = sort(df0['pv'].values)
        ok = _subsample(pv)

        qnull = -log10((0.5 + arange(len(pv))) / len(pv))
        qemp = -log10(pv)

        ax.plot(qnull[ok], qemp[ok], '-o', label=label, **style.get(label))

        if label not in style:
            style[label] = dict()

    ax.plot([0, qnull.max()], [0, qnull.max()], 'r')
    _plot_confidence_band(ok, qnull, alpha, ax)
    _set_axis_labels(ax)

    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')

    return ax


def _set_axis_labels(ax):
    ax.set_ylabel('-log$_{10}$pv observed')
    ax.set_xlabel('-log$_{10}$pv expected')
    ax.legend(loc='best')


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
    n = sum(nok)
    i = 0
    while i < len(idx) and idx[i] == n:
        i += 1
    idx = idx[i:]

    ok[where(nok)[0][idx]] = True

    ok[0] = True
    ok[-1] = True

    return ok


def _normalise_data(data):
    from pandas import DataFrame, concat

    if not isinstance(data, DataFrame):

        pvs = atleast_2d(data)

        dfs = []
        for i, pv in enumerate(pvs):
            df = DataFrame(columns=['pv', 'label'])
            df['pv'] = pv
            df['label'] = 'label{}'.format(i)
            dfs.append(df)

        return concat(dfs, ignore_index=True)

    return data
