from .._display import session_line as _session_line


def read(filepath, verbose=True):
    from numpy import load

    with _session_line("Reading {}...".format(filepath), disable=not verbose):
        return load(filepath)


def save(filepath, X, verbose=True):
    with _session_line("Saving {}...".format(filepath), disable=not verbose):
        save(filepath, X)


def _see(filepath, verbose=True):
    print(read(filepath, verbose=verbose))
