import os
import sys

PY2 = sys.version_info < (3, )

if PY2:
    from urllib import urlretrieve
    from urlparse import urlparse
else:
    from urllib.request import urlretrieve
    from urllib.parse import urlparse


def download(url, verbose=True):
    if verbose:
        print("Downloading {}...".format(url))
    urlretrieve(url, _filename(url))


def _filename(url):
    a = urlparse(url)
    return os.path.basename(a.path)
