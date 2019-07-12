def pytest_configure(*_):
    setup_tests_baseline()


def pytest_runtest_setup():
    setup_tests_baseline()


def setup_tests_baseline():
    _compatibility()
    import doctest
    import numpy as np
    import matplotlib as mpl

    mpl.use("agg")

    _pandas_format()
    doctest.ELLIPSIS_MARKER = "-ignore-"
    np.set_printoptions(precision=8, floatmode="fixed")


def _pandas_format():
    import pandas as pd

    pd.set_option("display.width", 88)
    pd.set_option("display.max_columns", 79)
    pd.set_option("display.max_rows", 60)
    pd.set_option("display.large_repr", "truncate")
    pd.set_option("display.float_format", "{:8.5f}".format)


def _compatibility():
    import warnings

    warnings.filterwarnings("ignore", message="numpy.dtype size changed")
    warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
