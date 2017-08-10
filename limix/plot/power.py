from __future__ import division

from numpy import argsort, asarray, linspace


def plot_power(df, style=None, ax=None):
    r"""Plot number of hits across significance levels.

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        Data frame with `pv` and `label` columns.
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

        from limix.plot import plot_power
        from pandas import DataFrame
        from numpy.random import RandomState
        from matplotlib import pyplot as plt
        random = RandomState(1)
        nsnps = 10000

        pv0 = list(random.rand(nsnps))
        pv1 = list(0.7 * random.rand(nsnps))

        fig = plt.figure(1, figsize=(5,5))
        plt.subplot(111)

        data = dict(pv=pv0 + pv1,
                    label=['label0'] * nsnps + ['label1'] * nsnps)
        df = DataFrame(data=data)
        plot_power(df)
        plt.show()
    """

    import matplotlib.pyplot as plt

    ax = plt.gca() if ax is None else ax

    labels = list(df['label'].unique())
    if style is None:
        style = {label: dict() for label in labels}

    alphas, nhits = _collect_nhits(df)

    for label in labels:
        ax.plot(
            alphas, asarray(nhits[label], int), label=label, **style.get(label))

    _set_labels(ax)

    return ax


def _collect_nhits(df):
    labels = list(df['label'].unique())
    alphas = linspace(0.01, 0.5, 500)
    nhits = {l: [] for l in labels}

    for i in range(len(alphas)):
        alpha = alphas[i]
        for label in labels:
            ix = df['label'] == label
            df_ = df.loc[ix, :]
            ntests = df_.shape[0]
            n = (df_['pv'] < alpha).sum()
            nhits[label] += [n]

    for label in labels:
        nhits[label] = asarray(nhits[label], int)

    return (alphas, nhits)


def _set_labels(ax):
    ax.set_xlabel('significance level')
    ax.set_ylabel('number of hits')
    ax.legend()


def plot_power_known(df, alpha=0.05, style=None, ax=None):

    from limix.stats import confusion_matrix
    import matplotlib.pyplot as plt

    ax = plt.gca() if ax is None else ax

    labels = list(df['label'].unique())
    if style is None:
        style = {label: dict() for label in labels}

    for label in labels:

        df0 = df.query("label=='%s'" % label)
        cm = confusion_matrix(df0)
        y = cm.tpr[1:]
        x = cm.fpr[1:]

        idx = argsort(x)
        x = x[idx]
        y = y[idx]

        ok = x <= alpha
        x = x[ok]
        y = y[ok]

        ax.plot(x, y * 100, label=label, **style.get(label))

    ax.set_xlabel('significance level')
    ax.set_ylabel('percentage of hits')
    ax.legend()

    return ax
