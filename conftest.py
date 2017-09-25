def get_pkg_name():
    from setuptools import find_packages

    return find_packages()[0]

collect_ignore = ["setup.py", "%s/_test.py" % get_pkg_name()]
