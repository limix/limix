import os
import pytest
import _pytest


def pytest_sessionstart(session):
    import pandas as pd

    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 79)
    pd.set_option("display.max_rows", 60)
    pd.set_option("display.large_repr", "truncate")


@pytest.fixture(autouse=True)
def _docdir(request):

    # Trigger ONLY for the doctests.
    doctest_plugin = request.config.pluginmanager.getplugin("doctest")
    if isinstance(request.node, doctest_plugin.DoctestItem):

        # Get the fixture dynamically by its name.
        tmpdir = request.getfuncargvalue("tmpdir")

        # Chdir only for the duration of the test.
        olddir = os.getcwd()
        tmpdir.chdir()
        yield
        os.chdir(olddir)

    else:
        # For normal tests, we have to yield, since this is a yield-fixture.
        yield
