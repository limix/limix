from __future__ import division

import warnings

from numpy import argsort, asarray, inf, isnan, percentile, clip


def plot_kinship(K, nclusters=1, style=None, ax=None):
    r"""Plot Kinship matrix.

    Parameters
    ----------
    K : array_like
        Kinship matrix.
    nclusters : int or str
        Number of blocks to be seen from the heatmap. It defaults to ``1``,
        which means that no ordering is performed. Pass 'auto' to automatically
        determine the number of clusters. Pass an integer to select the number
        of clusters.
    style : dict
        Keyword arguments forwarded to the :func:`matplotlib.axes.Axes.imshow`
        function.
    ax : :class:`matplotlib.axes.Axes`
        The target handle for this figure. If ``None``, the current axes is
        set.

    Returns
    -------
    :class:`matplotlib.axes.Axes`
        Axes.
    """

    import matplotlib.pyplot as plt

    ax = plt.gca() if ax is None else ax

    if style is None:
        style = dict()

    K = asarray(K, float)
    if nclusters == 'auto':
        K = _infer_clustering(K)
    elif nclusters > 1:
        K = _clustering(K, nclusters)

    cmin = percentile(K, 2)
    cmax = percentile(K, 98)
    K = clip(K, cmin, cmax)
    K = (K - K.min()) / (K.max() - K.min())

    mesh = ax.pcolormesh(K, cmap='RdBu_r')

    ax.set_aspect("equal")
    ax.set(xlim=(0, K.shape[1]), ylim=(0, K.shape[0]))
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    cb = ax.figure.colorbar(mesh, None, ax)

    #    ax.xaxis.set_ticks_position('both')
    #    ax.yaxis.set_ticks_position('both')
    #    ax.spines['right'].set_visible(True)
    #    ax.spines['top'].set_visible(True)
    #    ax.spines['left'].set_visible(True)
    #    ax.spines['bottom'].set_visible(True)
    #
    #    fig = plt.gcf()
    #    cb = fig.axes[-1]
    #    cb.xaxis.set_ticks_position('both')
    #    cb.yaxis.set_ticks_position('both')
    #    cb.outline.set_linewidth(1.0)
    #    #for spine in cb.spines.values():
    #    #    spine.set_visible(True)
    #    #    spine.set_linewidth(10.0)
    #    #    spine.set_color('black')
    #    #    spine.set_edgecolor('black')

    return ax


def _infer_clustering(K):
    from sklearn.metrics import silhouette_score

    scores = []
    nclusterss = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

    for nclusters in nclusterss:
        labels = _cluster(K, nclusters)
        # idx = argsort(labels)

        s = silhouette_score(K, labels, metric='correlation')
        scores.append(s)

    smallest = inf
    nclusters = -1
    for i in range(1, len(nclusterss)):
        d = scores[i] - scores[i - 1]
        if d < smallest:
            smallest = d
            nclusters = nclusterss[i - 1]

    return _clustering(K, nclusters)


def _cluster(K, n):
    from sklearn.cluster import SpectralClustering

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = SpectralClustering(n_clusters=n)
        m.fit(K)

    return m.labels_


def _clustering(K, n):
    idx = argsort(_cluster(K, n))
    return K[idx, :][:, idx]
