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
def download(url, verbose=True):
    filename = _filename(url)
    if os.path.exists(filename):
        if verbose:
            print("File {} already exists.".format(filename))
        return

    if verbose:
        print("Downloading {}...".format(url))
    urlretrieve(url, filename)


def _filename(url):
    a = urlparse(url)
    return os.path.basename(a.path)
