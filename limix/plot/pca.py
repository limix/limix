from __future__ import division

from numpy import asarray


def plot_pca(X, style=None, ax=None):
    r"""Plot the first two principal components of a design matrix.

    Parameters
    ----------
    X : array_like
        Design matrix.
    style : dict
        Keyword arguments forwarded to the :func:`matplotlib.pyplt.scatter`
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

        from numpy.random import RandomState
        import limix

        random = RandomState(0)
        X = random.randn(30, 10)

        limix.plot.plot_pca(X)
        limix.plot.show()
    """

    import matplotlib.pyplot as plt
    from sklearn import decomposition

    ax = plt.gca() if ax is None else ax

    if style is None:
        style = dict()

    X = asarray(X, float)

    pca = decomposition.PCA(n_components=2)
    pca.fit(X)
    X = pca.transform(X)

    ax.scatter(X[:, 0], X[:, 1])
    # ax.scatter(X[:, 0], X[:, 1], cmap=plt.cm.spectral)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    ax.set_xlabel('first principal component')
    ax.set_ylabel('second principal component')

    return ax
