import os
import pytest


def pytest_sessionstart(session):
    import doctest
    import matplotlib as mpl

    _compatibility()
    import pandas as pd

    mpl.use("agg")
    doctest.ELLIPSIS_MARKER = "-ignore-"

    pd.set_option("display.width", 88)
    pd.set_option("display.max_columns", 79)
    pd.set_option("display.max_rows", 60)
    pd.set_option("display.large_repr", "truncate")
    pd.set_option("display.float_format", "{:8.5f}".format)


@pytest.fixture(autouse=True)
def _docdir(request):

    # Trigger ONLY for the doctests or doctestplus.
    plug = request.config.pluginmanager.getplugin("doctest")
    if plug is None:
        plug = request.config.pluginmanager.getplugin("doctestplus")
        if plug is None:
            item = None
        else:
            item = plug._doctest_textfile_item_cls
    else:
        item = plug.DoctestItem

    if isinstance(request.node, item):
        # Get the fixture dynamically by its name.
        tmpdir = request.getfixturevalue("tmpdir")

        # Chdir only for the duration of the test.
        olddir = os.getcwd()
        tmpdir.chdir()
        yield
        os.chdir(olddir)

    else:
        # For normal tests, we have to yield, since this is a yield-fixture.
        yield


def _compatibility():
    import warnings

    warnings.filterwarnings("ignore", message="numpy.dtype size changed")
    warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
