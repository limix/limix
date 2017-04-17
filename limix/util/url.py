import sys

PY2 = sys.version_info < (3, )

if PY2:
    from urllib import urlretrieve
else:
    from urllib.request import urlretrieve
