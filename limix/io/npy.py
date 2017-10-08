from numpy import load


def see_kinship(filepath):
    # TODO: document
    import limix

    K = load(filepath)
    limix.plot.plot_kinship(K)
