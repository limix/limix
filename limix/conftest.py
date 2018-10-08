from __future__ import unicode_literals


def pytest_sessionstart(*args, **kwargs):
    import matplotlib as mpl

    mpl.use("agg")

    _compatibility()
    _pandas_format()


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
