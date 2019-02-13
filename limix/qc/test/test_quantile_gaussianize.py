import dask.array as da
import dask.dataframe as dd
import pandas as pd
import pytest
import xarray as xr
from limix.qc import quantile_gaussianize
from numpy import asarray, nan
from numpy.random import RandomState
from numpy.testing import assert_, assert_allclose


@pytest.fixture
def data2d():
    random = RandomState(0)
    X = random.randn(5, 2)

    R = asarray(
        [
            [0.430_727_299, 0.0],
            [0.0, 0.967_421_566],
            [0.967_421_566, -0.967_421_566],
            [-0.430_727_299, -0.430_727_299],
            [-0.967_421_566, 0.430_727_299],
        ]
    )
    Rt = asarray(
        [
            [0.430_727_299, -0.430_727_299],
            [-0.430_727_299, 0.430_727_299],
            [0.430_727_299, -0.430_727_299],
            [0.430_727_299, -0.430_727_299],
            [-0.430_727_299, 0.430_727_299],
        ]
    )

    Xnan = X.copy()
    Xnan[0, 0] = nan

    Rnan = asarray(
        [
            [nan, 0.0],
            [0.253_347_103, 0.967_421_566],
            [0.841_621_234, -0.967_421_566],
            [-0.253_347_103, -0.430_727_299],
            [-0.841_621_234, 0.430_727_299],
        ]
    )

    Rnant = asarray(
        [
            [nan, 0.0],
            [-0.430_727_299, 0.430_727_299],
            [0.430_727_299, -0.430_727_299],
            [0.430_727_299, -0.430_727_299],
            [-0.430_727_299, 0.430_727_299],
        ]
    )

    samples = [f"snp{i}" for i in range(X.shape[0])]
    candidates = [f"snp{i}" for i in range(X.shape[1])]

    return {
        "X": X,
        "Xnan": Xnan,
        "R": R,
        "Rnan": Rnan,
        "Rt": Rt,
        "Rnant": Rnant,
        "samples": samples,
        "candidates": candidates,
    }


@pytest.fixture
def data1d():
    X = asarray([0.0, 0.400_157_208_367_223_3, 0.978_737_984_105_739_2])
    Xnan = asarray([nan, 0.400_157_208_367_223_3, 0.978_737_984_105_739_2])
    R = asarray([-0.674_489_75, 0.0, 0.674_489_75])
    Rnan = asarray([nan, -0.430_727_299, 0.430_727_299])
    Rt = asarray([0.0, 0.0, 0.0])
    Rnant = [nan, 0, 0]
    samples = [f"snp{i}" for i in range(X.shape[0])]
    return {
        "X": X,
        "Xnan": Xnan,
        "R": R,
        "Rnan": Rnan,
        "Rt": Rt,
        "Rnant": Rnant,
        "samples": samples,
    }


def test_quantile_gaussianize_ndarray_1d(data1d):
    _assert_quantile_gaussianize(data1d, lambda X, _: asarray(X).copy())
    _assert_quantile_gaussianize_inplace(data1d, lambda X, _: asarray(X).copy())


def test_quantile_gaussianize_ndarray(data2d):
    _assert_quantile_gaussianize(data2d, lambda X, _: asarray(X).copy())
    _assert_quantile_gaussianize_inplace(data2d, lambda X, _: asarray(X).copy())


def test_quantile_gaussianize_pandas_series(data1d):
    _assert_quantile_gaussianize(data1d, lambda X, _: pd.Series(X.copy()))
    _assert_quantile_gaussianize_inplace(data1d, lambda X, _: pd.Series(X.copy()))


def test_quantile_gaussianize_pandas_dataframe(data2d):
    _assert_quantile_gaussianize(
        data2d, lambda X, d: pd.DataFrame(X.copy(), columns=d["candidates"])
    )
    _assert_quantile_gaussianize_inplace(
        data2d, lambda X, d: pd.DataFrame(X.copy(), columns=d["candidates"])
    )


def test_quantile_gaussianize_dask_array(data2d):
    _assert_quantile_gaussianize(data2d, lambda X, d: da.from_array(X.copy(), chunks=2))


def test_quantile_gaussianize_dask_dataframe(data2d):
    _assert_quantile_gaussianize(
        data2d, lambda X, d: dd.from_array(X.copy(), columns=d["candidates"])
    )


def test_quantile_gaussianize_xarray_dataarray(data2d):
    _assert_quantile_gaussianize(data2d, lambda X, _: xr.DataArray(X.copy()))
    _assert_quantile_gaussianize_inplace(data2d, lambda X, _: xr.DataArray(X.copy()))


def _assert_quantile_gaussianize(data, astype):
    qg = quantile_gaussianize

    X = astype(data["X"], data).copy()
    Xnan = astype(data["Xnan"], data).copy()
    assert_allclose(qg(X), data["R"])
    assert_allclose(X, astype(data["X"], data))
    assert_allclose(qg(Xnan), data["Rnan"])
    assert_allclose(Xnan, astype(data["Xnan"], data))

    assert_(isinstance(qg(X), type(astype(data["X"], data))))
    assert_(isinstance(qg(Xnan), type(astype(data["Xnan"], data))))

    X = astype(data["X"], data).copy()
    Xnan = astype(data["Xnan"], data).copy()
    assert_allclose(qg(X, axis=0), data["Rt"])
    assert_allclose(X, astype(data["X"], data))
    assert_allclose(qg(Xnan, axis=0), data["Rnant"])
    assert_allclose(Xnan, astype(data["Xnan"], data))

    assert_(isinstance(qg(X, axis=0), type(astype(data["X"], data))))
    assert_(isinstance(qg(Xnan, axis=0), type(astype(data["Xnan"], data))))


def _assert_quantile_gaussianize_inplace(data, astype):
    qg = quantile_gaussianize

    X = astype(data["X"], data).copy()
    Xnan = astype(data["Xnan"], data).copy()
    assert_allclose(qg(X, inplace=True), data["R"])
    assert_allclose(X, data["R"])
    assert_allclose(qg(Xnan, inplace=True), data["Rnan"])
    assert_allclose(Xnan, data["Rnan"])

    assert_(isinstance(qg(X, inplace=True), type(astype(data["X"], data))))
    assert_(isinstance(qg(Xnan, inplace=True), type(astype(data["Xnan"], data))))

    X = astype(data["X"], data).copy()
    Xnan = astype(data["Xnan"], data).copy()
    assert_allclose(qg(X, axis=0, inplace=True), data["Rt"])
    assert_allclose(X, data["Rt"])
    assert_allclose(qg(Xnan, axis=0, inplace=True), data["Rnant"])
    assert_allclose(Xnan, data["Rnant"])

    assert_(isinstance(qg(X, axis=0, inplace=True), type(astype(data["X"], data))))
    assert_(
        isinstance(qg(Xnan, axis=0, inplace=True), type(astype(data["Xnan"], data)))
    )
