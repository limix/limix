from .._display import session_line


def download(url, dest=None, verbose=True):
    """
    Download file.

    Parameters
    ----------
    url : str
        Url to the file.
    dest : str, optional
        File destination. The current working directory is used if ``None`` is passed.
        Defaults to ``None``.
    verbose : bool, optional
        ``True`` for displaying progress. Defaults to ``True``.

    Returns
    -------
    filepath : str
        File path to the downloaded file.
    """
    import os
    from urllib.request import urlretrieve

    if dest is None:
        dest = os.getcwd()

    filepath = os.path.join(dest, _filename(url))

    with session_line(f"Downloading {url}... ", disable=not verbose):
        urlretrieve(url, filepath)

    return filepath


def _filename(url):
    import os
    from urllib.parse import urlparse

    a = urlparse(url)
    return os.path.basename(a.path)
