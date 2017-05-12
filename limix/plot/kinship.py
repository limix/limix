from numpy import asarray


def plot_kinship(K):
    r"""Plot Kinship matrix.

    Parameters
    ----------
    K : array_like
        Kinship matrix.

    Examples
    --------

    .. plot::

        from limix.io.npy import see_kinship
        from limix.io.examples import numpy_kinship_file_example

        see_kinship(numpy_kinship_file_example())
    """
    from matplotlib import pyplot as plt
    from matplotlib import cm as cm

    K = asarray(K, float)
    mi = K.min()
    ma = K.max()
    K = (K - mi) / (ma - mi)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    cax = ax1.imshow(K, interpolation='nearest')
    fig.colorbar(cax)
    plt.show()
