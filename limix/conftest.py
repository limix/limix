from __future__ import unicode_literals


def pytest_sessionstart(*args, **kwargs):
    import matplotlib as mpl

    mpl.use("agg")

    _compatibility()
    pandas_format()


@pytest.fixture(autouse=True)
def _docdir(request):

    # Trigger ONLY for the doctests.
    plug = request.config.pluginmanager.getplugin("doctest")
    if plug is not None and isinstance(request.node, plug.DoctestItem):

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
