def get_pkg_name():
    from setuptools import find_packages

    return find_packages()[0]


collect_ignore = [
    "doc/conf.py", "setup.py", "{}/testit.py".format(get_pkg_name())
]
