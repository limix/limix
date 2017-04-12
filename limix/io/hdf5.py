import h5py
import dask.array as da
import dask.dataframe as dd
import pandas as pd

from numpy.testing import assert_allclose, assert_equal

class h5data_fetcher(object):
    def __init__(self, filename):
        self._filename = filename

    def __enter__(self):
        self._f = h5py.File(self._filename, 'r')
        return self

    def fetch(self, data_path):
        data = self._f[data_path]
        if data.chunks is None:
            chunks = data.shape
        else:
            chunks = data.chunks
        return da.from_array(data, chunks=chunks)

    def __exit__(self, *exc):
        self._f.close()
