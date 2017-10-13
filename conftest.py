import matplotlib
matplotlib.use('Agg')

import doctest
from limix import npy_doctest

doctest.OutputChecker = npy_doctest.FlexNumOutputChecker


def get_pkg_name():
    from setuptools import find_packages
    return find_packages()[0]


collect_ignore = [
    "conftest.py", "doc/conf.py", "setup.py",
    "{}/testit.py".format(get_pkg_name())
]
