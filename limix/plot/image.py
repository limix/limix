def see_image(filepath):
    import matplotlib.pyplot as plt

    plt.imshow(plt.imread(filepath))
    axes = plt.gca()
    plt.tight_layout()
    axes.set_position([0, 0, 1, 1])
    axes.xaxis.set_visible(False)
    axes.yaxis.set_visible(False)
    plt.show()
