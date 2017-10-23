from limix.util import Timer


def see_image(filepath, verbose):
    r"""Refer to :method:`limix.LimixPlot.see_image`."""
    import matplotlib.pyplot as plt

    with Timer(desc="Reading %s..." % filepath, disable=not verbose):
        plt.imshow(plt.imread(filepath))
        axes = plt.gca()
        plt.tight_layout()
        axes.set_position([0, 0, 1, 1])
        axes.xaxis.set_visible(False)
        axes.yaxis.set_visible(False)
    plt.show()
