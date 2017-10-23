import os
import sys

PY2 = sys.version_info < (3, )

if PY2:
    from urllib import urlretrieve
    from urlparse import urlparse
else:
    from urllib.request import urlretrieve
    from urllib.parse import urlparse


# TODO: document it
def download(url, dest=None, verbose=True):
    if dest is None:
        dest = os.getcwd()

    filepath = os.path.join(dest, _filename(url))
    if os.path.exists(filepath):
        if verbose:
            print("File {} already exists.".format(filepath))
        return

    if verbose:
        print("Downloading {}...".format(url))
    urlretrieve(url, filepath)


def _filename(url):
    a = urlparse(url)
    return os.path.basename(a.path)
