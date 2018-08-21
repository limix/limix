from ..display import timer_text


def see(filepath, verbose=True):
    # TODO: document
    print(read(filepath, verbose=verbose))


def read(filepath, verbose=True):
    from numpy import load

    with timer_text("Reading {}...".format(filepath), disable=not verbose):
        return load(filepath)


def save(filepath, X, verbose=True):
    with timer_text("Saving {}...".format(filepath), disable=not verbose):
        save(filepath, X)
