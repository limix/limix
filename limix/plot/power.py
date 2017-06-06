from __future__ import division

from numpy import asarray, linspace


def plot_power_curve(df, color=None, ax=None):
    r"""Plot number of hits across significance levels.

    Parameters
    ----------

    df : :class:`pandas.DataFrame`
        Data frame with `pv` and `label` columns.
    color : dict
        Map colors to labels.
    ax : :class:`matplotlib.axes.AxesSubplot`
        The target handle for this figure. If None, the current axes is set.

    Returns
    -------
    :class:`matplotlib.axes.AxesSubplot`
        Axes.

    Examples
    --------

    .. plot::

        from limix.plot import plot_power_curve
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
        plot_power_curve(df)
        plt.show()
    """

    import matplotlib.pyplot as plt

    ax = plt.gca() if ax is None else ax
    labels = list(df['label'].unique())

    if color is None:
        colors = _get_default_colors()
        color = {m: colors[i] for (i, m) in enumerate(labels)}

    alphas, nhits = _collect_nhits(df)

    for label in labels:
        ax.plot(
            alphas,
            asarray(nhits[label], int),
            color=color[label],
            label=label)

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


def _get_default_colors():
    return ['red', 'green', 'blue']
