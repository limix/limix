import dask.dataframe as dd
import pandas as pd
import pytest
import xarray as xr
from limix.qc import quantile_gaussianize
from numpy import asarray, nan, ndarray
from numpy.random import RandomState
from numpy.testing import assert_, assert_allclose


@pytest.fixture
def data():
    random = RandomState(0)
    X0 = random.randint(0, 3, size=(5, 2))
    X1 = random.randn(5, 2)

    R0 = asarray(
        [
            [-0.210_428_394_247_924_84, -0.210_428_394_247_924_84],
            [-0.210_428_394_247_924_84, -0.210_428_394_247_924_84],
            [0.967_421_566_101_701, 0.674_489_750_196_081_7],
            [-0.210_428_394_247_924_84, 0.674_489_750_196_081_7],
            [-0.210_428_394_247_924_84, -0.967_421_566_101_701],
        ]
    )
    R1 = asarray(
        [
            [-0.967_421_566_101_701, -0.967_421_566_101_701],
            [-0.430_727_299_295_457_44, 0.967_421_566_101_701],
            [0.430_727_299_295_457_56, 0.0],
            [0.967_421_566_101_701, -0.430_727_299_295_457_44],
            [0.0, 0.430_727_299_295_457_56],
        ]
    )
    X1nan = X1.copy()
    X1nan[0, 0] = nan
    R1nan = asarray(
        [
            [nan, -0.967_421_566_101_701],
            [-0.841_621_233_572_914_3, 0.967_421_566_101_701],
            [0.253_347_103_135_799_7, 0.0],
            [0.841_621_233_572_914_2, -0.430_727_299_295_457_44],
            [-0.253_347_103_135_799_7, 0.430_727_299_295_457_56],
        ]
    )

    samples = [f"snp{i}" for i in range(X0.shape[0])]
    candidates = [f"snp{i}" for i in range(X0.shape[1])]

    return {
        "X0": X0,
        "X1": X1,
        "X1nan": X1nan,
        "R0": R0,
        "R1": R1,
        "R1nan": R1nan,
        "samples": samples,
        "candidates": candidates,
    }


def test_quantile_gaussianize_ndarray(data):
    assert_allclose(quantile_gaussianize(data["X0"]), data["R0"])
    assert_allclose(quantile_gaussianize(data["X1"]), data["R1"])
    assert_allclose(quantile_gaussianize(data["X1nan"]), data["R1nan"])
    assert_(isinstance(quantile_gaussianize(data["X0"]), ndarray))
    assert_(isinstance(quantile_gaussianize(data["X1"]), ndarray))
    assert_(isinstance(quantile_gaussianize(data["X1nan"]), ndarray))


def test_quantile_gaussianize_pandas_dataframe(data):
    def get(X):
        return pd.DataFrame(X, columns=data["candidates"])

    assert_allclose(quantile_gaussianize(get(data["X0"])), data["R0"])
    assert_allclose(quantile_gaussianize(get(data["X1"])), data["R1"])
    assert_allclose(quantile_gaussianize(get(data["X1nan"])), data["R1nan"])
    assert_(isinstance(quantile_gaussianize(get(data["X0"])), pd.DataFrame))
    assert_(isinstance(quantile_gaussianize(get(data["X1"])), pd.DataFrame))
    assert_(isinstance(quantile_gaussianize(get(data["X1nan"])), pd.DataFrame))


def test_quantile_gaussianize_dask_dataframe(data):
    def get(X):
        df = pd.DataFrame(X, columns=data["candidates"])
        return dd.from_pandas(df, npartitions=2)

    assert_allclose(quantile_gaussianize(get(data["X0"])), data["R0"])
    assert_allclose(quantile_gaussianize(get(data["X1"])), data["R1"])
    assert_allclose(quantile_gaussianize(get(data["X1nan"])), data["R1nan"])
    assert_(isinstance(quantile_gaussianize(get(data["X0"])), dd.DataFrame))
    assert_(isinstance(quantile_gaussianize(get(data["X1"])), dd.DataFrame))
    assert_(isinstance(quantile_gaussianize(get(data["X1nan"])), dd.DataFrame))


def test_quantile_gaussianize_xarray_dataarray(data):
    def get(X):
        return xr.DataArray(
            X,
            coords={"candidates": data["candidates"], "samples": data["samples"]},
            dims=["samples", "candidates"],
        )

    assert_allclose(quantile_gaussianize(get(data["X0"])), data["R0"])
    assert_allclose(quantile_gaussianize(get(data["X1"])), data["R1"])
    assert_allclose(quantile_gaussianize(get(data["X1nan"])), data["R1nan"])
    assert_(isinstance(quantile_gaussianize(get(data["X0"])), xr.DataArray))
    assert_(isinstance(quantile_gaussianize(get(data["X1"])), xr.DataArray))
    assert_(isinstance(quantile_gaussianize(get(data["X1nan"])), xr.DataArray))
