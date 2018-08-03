import os
import sys
import warnings

PY2 = sys.version_info < (3,)

if PY2:
    from urllib import urlretrieve
    from urlparse import urlparse
else:
    from urllib.request import urlretrieve
    from urllib.parse import urlparse


# TODO: document it
def download(url, dest=None, verbose=True, force=False):
    if dest is None:
        dest = os.getcwd()

    filepath = os.path.join(dest, _filename(url))
    if not force and os.path.exists(filepath):
        warnings.warn("File {} already exists.".format(filepath))
        print("Set `force` to `True` in order to overwrite the existing file.")
        return

    if verbose:
        print("Downloading {}...".format(url))
    urlretrieve(url, filepath)


def _filename(url):
    a = urlparse(url)
    return os.path.basename(a.path)
