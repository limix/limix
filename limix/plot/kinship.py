from numpy import asarray


def plot_kinship(K, ax=None):
    r"""Plot Kinship matrix.

    Parameters
    ----------
    K : array_like
        Kinship matrix.
    ax : :class:`matplotlib.axes.AxesSubplot`
        The target handle for this figure. If None, the current axes is set.

    Returns
    -------
    :class:`matplotlib.axes.AxesSubplot`
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

    K = asarray(K, float)
    mi = K.min()
    ma = K.max()
    K = (K - mi) / (ma - mi)

    cax = ax.imshow(K, interpolation='nearest')
    ax.figure.colorbar(cax)

    return ax
