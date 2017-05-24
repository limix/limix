from numpy import load


def see_kinship(filepath):
    import limix

    K = load(filepath)
    limix.plot.plot_kinship(K)
