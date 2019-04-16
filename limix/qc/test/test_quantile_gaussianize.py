import pytest
from numpy import asarray, nan
from numpy.random import RandomState

from limix.qc import quantile_gaussianize as qg

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
    random = RandomState(0)
    X = random.randn(5, 2)
    X[0, 0] = nan

    R = asarray(
        [
            [nan, 0.0],
            [0.253_347_103, 0.967_421_566],
            [0.841_621_234, -0.967_421_566],
            [-0.253_347_103, -0.430_727_299],
            [-0.841_621_234, 0.430_727_299],
        ]
    )

    Rt = asarray(
        [
            [nan, 0.0],
            [-0.430_727_299, 0.430_727_299],
            [0.430_727_299, -0.430_727_299],
            [0.430_727_299, -0.430_727_299],
            [-0.430_727_299, 0.430_727_299],
        ]
    )

    samples = [f"sample{i}" for i in range(X.shape[0])]
    candidates = [f"snp{i}" for i in range(X.shape[1])]

    return {"X": X, "R": R, "Rt": Rt, "samples": samples, "candidates": candidates}


@pytest.fixture
def data1d():
    X = asarray([nan, 0.400_157_208_367_223_3, 0.978_737_984_105_739_2])
    R = asarray([nan, -0.430_727_299, 0.430_727_299])
    Rt = [nan, 0, 0]
    samples = [f"sample{i}" for i in range(X.shape[0])]
    return {"X": X, "R": R, "Rt": Rt, "samples": samples, "candidates": None}


def test_impute_ndarray_1d(data1d):
    assert_ndarray_1d(qg, data1d)


def test_impute_ndarray_2d(data2d):
    assert_ndarray_2d(qg, data2d)


def test_impute_pandas_series(data1d):
    assert_pandas_series(qg, data1d)


def test_impute_pandas_dataframe(data2d):
    assert_pandas_dataframe(qg, data2d)


def test_impute_dask_array(data2d):
    assert_dask_array(qg, data2d, test_inplace=False)


def test_impute_dask_dataframe(data2d):
    assert_dask_dataframe(qg, data2d)


def test_impute_xarray_dataarray(data2d):
    assert_xarray_dataarray(qg, data2d)
