from setuptools import find_packages

_pkg = find_packages()[0]

collect_ignore = ["setup.py", "%s/_test.py" % _pkg]
