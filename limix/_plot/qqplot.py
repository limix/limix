from __future__ import division

from numpy import (arange, ascontiguousarray, atleast_2d, flipud, inf, insert,
                   linspace, log10, ones, percentile, searchsorted, sort)
from numpy import sum as npsum
from numpy import where
from scipy.special import betaincinv


def qqplot(df, alpha, cutoff=0.1, style=None, ax=None, limits=None):
    r"""Quantile-Quantile plot of observed p-values versus theoretical ones.

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        Data frame.
    alpha : float
        Significance level defining the band boundary. Defaults to ``0.05``.
    cutoff : float
        P-values higher than `cutoff` will not be plotted.
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

    qmin = +inf
    qmax = -inf
    for label, df0 in df.groupby('label', sort=True):
        pv = sort(df0['pv'].values)
        ok = _subsample(pv, cutoff)

        qnull = -log10((0.5 + arange(len(pv))) / len(pv))

        qemp = -log10(pv)

        sty = style.get(label, {})
        ax.plot(
            qnull[ok], qemp[ok], 'o', markeredgecolor=None, label=label, **sty)

        qmin = min(qmin, min(qnull[ok].min(), qemp[ok].min()))
        qmax = max(qmax, max(qnull[ok].max(), qemp[ok].max()))

        if label not in style:
            style[label] = dict()

    if limits is not None:
        qmin = limits[0]
        qmax = limits[1]

    ax.plot([qmin, qmax], [qmin, qmax], color='black', zorder=0)
    _plot_confidence_band(ok, qnull, alpha, ax, qmax)
    _set_axis_labels(ax)

    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')

    ax.set_xlim(qmin, qmax)
    ax.set_ylim(qmin, qmax)

    ax.set_aspect('equal', 'box')
    ax.apply_aspect()

    ticks = ax.get_yticks()
    ymax = ax.get_xbound()[1]
    ticks = [t for t in ticks if t < ymax]

    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

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


def _plot_confidence_band(ok, null_qvals, significance_level, ax, qmax):

    (bo, to) = _rank_confidence_band(len(null_qvals), significance_level, ok)

    bo = -log10(bo)
    to = -log10(to)

    m = null_qvals[ok]
    m = insert(m, 0, qmax)

    d = m[0] - m[1]

    bo = insert(bo, 0, bo[0] + d * (bo[0] - bo[1]))
    to = insert(to, 0, to[0] + d * (to[0] - to[1]))

    ax.fill_between(
        m, bo, to, facecolor='#DDDDDD', linewidth=0, zorder=-1, alpha=1.0)


def _subsample(pvalues, cutoff):
    resolution = 500

    if len(pvalues) <= resolution:
        return ones(len(pvalues), dtype=bool)

    ok = pvalues <= percentile(pvalues, cutoff)
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
