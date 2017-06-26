from numpy import asarray


def plot_kinship(K, style=None, ax=None):
    r"""Plot Kinship matrix.

    Parameters
    ----------
    K : array_like
        Kinship matrix.
    style : dict
        Keyword arguments forwarded to the :func:`matplotlib.axes.Axes.imshow`
        function.
    ax : :class:`matplotlib.axes.Axes`
        The target handle for this figure. If None, the current axes is set.

    Returns
    -------
    :class:`matplotlib.axes.Axes`
        Axes.

    Examples
    --------

    .. plot::

        import numpy as np
        from matplotlib import pyplot as plt
        from limix.io.examples import numpy_kinship_file_example
        from limix.plot import plot_kinship

        K = np.load(numpy_kinship_file_example())
        plot_kinship(K)
        plt.tight_layout()
        plt.show()
    """
    import matplotlib.pyplot as plt

    ax = plt.gca() if ax is None else ax

    if style is None:
        style = dict()

    if 'interpolation' not in style:
        style['interpolation'] = 'nearest'

    K = asarray(K, float)
    mi = K.min()
    ma = K.max()
    K = (K - mi) / (ma - mi)

    cax = ax.imshow(K, **style)
    ax.figure.colorbar(cax)

    return ax
