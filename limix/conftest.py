from __future__ import unicode_literals

import pytest


def pytest_configure(*_):
    import doctest
    import matplotlib as mpl

    mpl.use("agg")

    _compatibility()
    pandas_format()
    doctest.ELLIPSIS_MARKER = "-ignore-"


@pytest.fixture(autouse=True)
def _docdir(request):
    import os

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


def pandas_format():
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
    warnings.filterwarnings(
        "ignore",
        message=(
            "Using or importing the ABCs from 'collections' instead"
            " of from 'collections.abc' is deprecated, and in 3.8 it will stop working"
        ),
    )
