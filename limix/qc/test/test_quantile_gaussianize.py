import dask.array as da
import dask.dataframe as dd
import pandas as pd
import pytest
import xarray as xr
from limix.qc import quantile_gaussianize as qg
from numpy import asarray, nan
from numpy.random import RandomState
from .util import assert_mat_proc, assert_mat_proc_inplace


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


def test_quantile_gaussianize_ndarray_1d(data1d):
    assert_mat_proc(qg, data1d, lambda X, *_: asarray(X).copy())
    assert_mat_proc_inplace(qg, data1d, lambda X, *_: asarray(X).copy())


def test_quantile_gaussianize_ndarray_2d(data2d):
    assert_mat_proc(qg, data2d, lambda X, *_: asarray(X).copy())
    assert_mat_proc_inplace(qg, data2d, lambda X, *_: asarray(X).copy())


def test_quantile_gaussianize_pandas_series(data1d):
    def get(X, samples, _):
        return pd.Series(X.copy(), index=samples)

    assert_mat_proc(qg, data1d, get)
    assert_mat_proc_inplace(qg, data1d, get)


def test_quantile_gaussianize_pandas_dataframe(data2d):
    def get(X, samples, candidates):
        return pd.DataFrame(X.copy(), index=samples, columns=candidates)

    assert_mat_proc(qg, data2d, get)
    assert_mat_proc_inplace(qg, data2d, get)


def test_quantile_gaussianize_dask_array(data2d):
    assert_mat_proc(qg, data2d, lambda X, *_: da.from_array(X.copy(), chunks=2))


def test_quantile_gaussianize_dask_dataframe(data2d):
    def get(X, samples, candidates):
        df = pd.DataFrame(X.copy(), index=samples, columns=candidates)
        df = dd.from_pandas(df, npartitions=2)
        return df

    assert_mat_proc(qg, data2d, get)


def test_quantile_gaussianize_xarray_dataarray(data2d):
    assert_mat_proc(qg, data2d, lambda X, *_: xr.DataArray(X.copy()))
    assert_mat_proc_inplace(qg, data2d, lambda X, *_: xr.DataArray(X.copy()))
