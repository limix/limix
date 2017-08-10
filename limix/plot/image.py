from limix.util import Timer


def see_image(filepath):
    r"""Plot a image represented in a file.

    It uses :func:`matplotlib.pyplot.imread` in order to read the provided
    file and plot.

    Parameters
    ----------
    filepath : str
        File path to the image.
    """
    import matplotlib.pyplot as plt

    with Timer(desc="Reading %s..." % filepath):
        plt.imshow(plt.imread(filepath))
        axes = plt.gca()
        plt.tight_layout()
        axes.set_position([0, 0, 1, 1])
        axes.xaxis.set_visible(False)
        axes.yaxis.set_visible(False)
    plt.show()
