from __future__ import division

from numpy import asarray


def plot_pca(X, style=None, ax=None):
    import matplotlib.pyplot as plt
    from sklearn import decomposition

    ax = plt.gca() if ax is None else ax

    if style is None:
        style = dict()

    X = asarray(X, float)

    pca = decomposition.PCA(n_components=2)
    pca.fit(X)
    X = pca.transform(X)

    ax.plot(X[:, 0], X[:, 1], 'o', markersize=4, **style)

    ax.set_xlabel('first component')
    ax.set_ylabel('second component')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')

    return ax
