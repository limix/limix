def get_pkg_name():
    from setuptools import find_packages

    return find_packages()[0]


import matplotlib

matplotlib.use('Agg')

collect_ignore = [
    "doc/conf.py", "setup.py", "{}/testit.py".format(get_pkg_name())
]
