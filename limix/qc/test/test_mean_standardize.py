import pytest
from numpy import array, nan
from numpy.random import RandomState

from limix.qc import mean_standardize

from .util import (
    assert_dask_array,
    assert_dask_dataframe,
    assert_ndarray_1d,
    assert_ndarray_2d,
    assert_pandas_dataframe,
    assert_pandas_series,
    assert_xarray_dataarray,
)


@pytest.fixture
def data2d():
    X = RandomState(0).randn(3, 2)
    X[0, 0] = nan
    R = array(
        [[nan, -0.117_142_199_8], [-1.0, 1.279_107_161_9], [1.0, -1.161_964_962_2]]
    )
    Rt = array([[nan, 0.0], [-1.0, 1.0], [1.0, -1.0]])

    samples = [f"sample{i}" for i in range(X.shape[0])]
    candidates = [f"snp{i}" for i in range(X.shape[1])]

    return {"X": X, "R": R, "Rt": Rt, "samples": samples, "candidates": candidates}


@pytest.fixture
def data1d():
    X = RandomState(0).randn(3)
    X[0] = nan
    R = array([nan, -1.0, 1.0])
    Rt = array([nan, 0.0, 0.0])

    samples = [f"sample{i}" for i in range(X.shape[0])]
    return {"X": X, "R": R, "Rt": Rt, "samples": samples, "candidates": None}


def test_mean_standardize_ndarray_1d(data1d):
    assert_ndarray_1d(mean_standardize, data1d)


def test_mean_standardize_ndarray_2d(data2d):
    assert_ndarray_2d(mean_standardize, data2d)


def test_mean_standardize_pandas_series(data1d):
    assert_pandas_series(mean_standardize, data1d)


def test_mean_standardize_pandas_dataframe(data2d):
    assert_pandas_dataframe(mean_standardize, data2d)


def test_mean_standardize_dask_array(data2d):
    assert_dask_array(mean_standardize, data2d, test_inplace=False)


def test_mean_standardize_dask_dataframe(data2d):
    assert_dask_dataframe(mean_standardize, data2d)


def test_mean_standardize_xarray_dataarray(data2d):
    assert_xarray_dataarray(mean_standardize, data2d)
