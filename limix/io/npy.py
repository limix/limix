from .._display import session_line as _session_line


def read(filepath, verbose=True):
    """
    Read NumPy arrays saved in a file.

    Parameters
    ----------
    filepath : str
        File path.
    verbose : bool, optional
        Defaults to ``True``.

    Returns
    -------
    array : array, tuple, dict, etc.
        Data stored in the file.
    """
    from numpy import load

    with _session_line("Reading {}...".format(filepath), disable=not verbose):
        return load(filepath)


def _see(filepath, verbose=True):
    print(read(filepath, verbose=verbose))
