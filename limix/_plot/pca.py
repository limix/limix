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

    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(X.min(), X.max())
    ax.set_aspect('equal', 'box')
    ax.apply_aspect()

    xbound = ax.get_xbound()
    xticks = ax.get_xticks()
    ax.set_yticks(xticks)
    ax.set_ybound(xbound)

    vxticks = [t for t in xticks if t >= xbound[0] and t <= xbound[1]]
    a = vxticks[0] - xbound[0]
    b = xbound[1] - vxticks[-1]
    if a > b:
        ax.set_xbound((xbound[0], vxticks[-1] + a))
    else:
        ax.set_xbound((vxticks[0] - n, xbound[1]))
    ax.set_ybound(ax.get_xbound())

    #ax.set_xlabel('first component')
    #ax.set_ylabel('second component')

    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')

    ax.grid(True, which='major', axis='both', alpha=1.0)

    return ax
