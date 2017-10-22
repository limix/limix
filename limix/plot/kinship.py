from __future__ import division

import warnings

from numpy import argsort, asarray, inf


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
    from seaborn import heatmap

    ax = plt.gca() if ax is None else ax

    if style is None:
        style = dict()

    K = asarray(K, float)
    if nclusters == 'auto':
        K = _infer_clustering(K)
    elif nclusters > 1:
        K = _clustering(K, nclusters)

    heatmap(
        K,
        ax=ax,
        linewidths=0,
        xticklabels=False,
        cmap='RdBu_r',
        yticklabels=False,
        square=True,
        robust=True,
        **style)

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
