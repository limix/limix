from numpy import load, save


def see_kinship(filepath):
    # TODO: document
    import limix

    K = load(filepath)
    limix.plot.plot_kinship(K)


def save_kinship(filepath, K, verbose=True):
    if verbose:
        print("Saving {}...".format(filepath))
    save(filepath, K)
