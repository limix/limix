from __future__ import division

from numpy import arange, asarray, cumsum, flipud, log10, arange, linspace


def plot_bonferroni(df, color=None, ax=None):

    import matplotlib.pyplot as plt

    ax = plt.gca() if ax is None else ax
    labels = list(df['label'].unique())

    if color is None:
        colors = _get_default_colors()
        color = {m:colors[i] for (i, m) in enumerate(labels)}

    alphas = linspace(0, 0.3, 30)
    y = {l:[] for l in labels}

    for i in range(len(alphas)):
        alpha = alphas[i]
        for label in labels:
            ix = df['label'] == label
            df_ = df.loc[ix, :]
            ntests = df_.shape[0]
            nhits = (df_['pv'] < alpha / ntests).sum()
            y[label] += [nhits]

    for label in labels:
        ax.plot(alphas, y[label], color=color[label], label=label)

    ax.set_xlabel('significance level')
    ax.set_ylabel('number of hits')
    ax.legend()
    return ax


# def plot_fdr(df, fdr=0.01, color=None, ax=None):
#
#     import matplotlib.pyplot as plt
#
#     ax = plt.gca() if ax is None else ax
#     labels = list(df['label'].unique())
#
#     if color is None:
#         colors = _get_default_colors()
#         color = {m:colors[i] for (i, m) in enumerate(labels)}
#
#     left = 0
#     for label in labels:
#         ix = df['label'] == label
#         df_ = df.loc[ix, :]
#         nd = (df_['pv'] < fdr).sum()
#         ax.bar(left, nd, color=color[label])
#         left += 1.0
#
#     ax.set_xticks(arange(0, left))
#     ax.set_xticklabels(labels)
#     ax.set_ylabel('number of hits (FDR %.1f%%)' % (100 * fdr))
#     return ax

def _get_default_colors():
    return ['red', 'green', 'blue']

if __name__ == '__main__':
    import numpy as np
    import pandas as pd
    labels = ['QEP'] * 100 + ['LMM'] * 100
    pvs = list(np.random.randn(100)) + list(np.random.randn(100))

    df = pd.DataFrame(data=dict(label=labels, pv=pvs))

    import matplotlib.pyplot as plt
    ax = plot_fdr(df)
    plt.show()
