import sys

from .._display import session_line

PY2 = sys.version_info < (3,)


# TODO: document it
def download(url, dest=None, verbose=True):
    import os

    if PY2:
        from urllib import urlretrieve
    else:
        from urllib.request import urlretrieve

    if dest is None:
        dest = os.getcwd()

    filepath = os.path.join(dest, _filename(url))

    with session_line("Downloading {}... ".format(url), disable=not verbose):
        urlretrieve(url, filepath)

    return filepath


def _filename(url):
    import os

    if PY2:
        from urlparse import urlparse
    else:
        from urllib.parse import urlparse

    a = urlparse(url)
    return os.path.basename(a.path)
