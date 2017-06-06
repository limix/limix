from __future__ import division

from numpy import asarray, linspace


def plot_power_curve(df, color=None, ax=None):

    import matplotlib.pyplot as plt

    ax = plt.gca() if ax is None else ax
    labels = list(df['label'].unique())

    if color is None:
        colors = _get_default_colors()
        color = {m: colors[i] for (i, m) in enumerate(labels)}

    alphas, nhits = _collect_nhits(df)

    for label in labels:
        ax.plot(
            alphas, asarray(y[label], int), color=color[label], label=label)

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

    nhits = asarray(nhits[label], int)

    return (alphas, nhits)


def _set_labels(ax):
    ax.set_xlabel('significance level')
    ax.set_ylabel('number of hits')
    ax.legend()


def _get_default_colors():
    return ['red', 'green', 'blue']
