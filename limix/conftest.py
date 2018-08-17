from __future__ import unicode_literals

import pkg_resources
import pytest
import os


def pytest_sessionstart(*args, **kwargs):
    pandas_format()


def pandas_format():
    import pandas as pd

    pd.set_option("display.width", 88)
    pd.set_option("display.max_columns", 79)
    pd.set_option("display.max_rows", 60)
    pd.set_option("display.large_repr", "truncate")
    pd.set_option("display.float_format", "{:8.5f}".format)


@pytest.fixture
def datadir(tmpdir, request):
    class DataDir(object):
        def __init__(self, tmpdir, module_path):
            self._tmpdir = tmpdir
            self._module_path = module_path

        def add(self, resource_path):
            resource_package = self._module_path
            name = resource_path.split("/")[-1]
            content = pkg_resources.resource_string(resource_package, resource_path)
            with open(os.path.join(self._tmpdir, name), "wb") as f:
                f.write(content)

        @property
        def tmpdir(self):
            return self._tmpdir

    return DataDir(tmpdir, request.module.__name__)
