from copy import copy

import warnings
import dask.array as da
import dask.dataframe as dd
import pandas as pd
import pytest
import xarray as xr
from limix.qc import mean_impute
from numpy import asarray, nan
from numpy.random import RandomState
from numpy.testing import assert_allclose


@pytest.fixture
def data2d():
    random = RandomState(0)
    X = random.randn(3, 4)
    X[0, 0] = nan

    R = asarray(
        [
            [
                0.882_169_569_178_204_8,
                0.400_157_208_367_223_3,
                0.978_737_984_105_739_2,
                2.240_893_199_201_458,
            ],
            [
                1.867_557_990_149_967_5,
                -0.977_277_879_876_411,
                0.950_088_417_525_589_4,
                -0.151_357_208_297_697_9,
            ],
            [
                -0.103_218_851_793_557_84,
                0.410_598_501_938_372_33,
                0.144_043_571_160_878,
                1.454_273_506_962_975,
            ],
        ]
    )
    Rt = asarray(
        [
            [
                1.206_596_130_558_14,
                0.400_157_208_367_223_3,
                0.978_737_984_105_739_2,
                2.240_893_199_201_458,
            ],
            [
                1.867_557_990_149_967_5,
                -0.977_277_879_876_411,
                0.950_088_417_525_589_4,
                -0.151_357_208_297_697_9,
            ],
            [
                -0.103_218_851_793_557_84,
                0.410_598_501_938_372_33,
                0.144_043_571_160_878,
                1.454_273_506_962_975,
            ],
        ]
    )

    samples = [f"snp{i}" for i in range(X.shape[0])]
    candidates = [f"snp{i}" for i in range(X.shape[1])]

    return {"X": X, "R": R, "Rt": Rt, "samples": samples, "candidates": candidates}


@pytest.fixture
def data1d():
    X = asarray([nan, 0.400_157_208_367_223_3, 0.978_737_984_105_739_2])
    R = asarray(
        [0.689_447_596_236_481_2, 0.400_157_208_367_223_3, 0.978_737_984_105_739_2]
    )
    Rt = R.copy()
    Rt[0] = nan
    samples = [f"snp{i}" for i in range(X.shape[0])]
    return {"X": X, "R": R, "Rt": Rt, "samples": samples}


def test_impute_ndarray_1d(data1d):
    _assert_impute(lambda d: d["X"].copy(), data1d)
    _assert_impute_inplace(lambda d: d["X"].copy(), data1d)


def test_impute_ndarray_2d(data2d):
    _assert_impute(lambda d: d["X"].copy(), data2d)
    _assert_impute_inplace(lambda d: d["X"].copy(), data2d)


def test_impute_pandas_series(data1d):
    _assert_impute(lambda d: pd.Series(d["X"]), data1d)
    _assert_impute_inplace(lambda d: pd.Series(d["X"]), data1d)


def test_impute_pandas_dataframe(data2d):
    _assert_impute(lambda d: pd.DataFrame(d["X"], columns=d["candidates"]), data2d)
    _assert_impute_inplace(
        lambda d: pd.DataFrame(d["X"], columns=d["candidates"]), data2d
    )


def test_impute_dask_array(data2d):
    _assert_impute(lambda d: da.from_array(d["X"].copy(), chunks=2), data2d)
    _assert_impute_inplace(lambda d: da.from_array(d["X"].copy(), chunks=2), data2d)


def test_impute_dask_dataframe(data2d):
    _assert_impute(lambda d: dd.from_array(d["X"].copy()), data2d)


def test_impute_xarray_dataarray(data2d):
    _assert_impute(lambda d: xr.DataArray(d["X"].copy()), data2d)
    _assert_impute_inplace(lambda d: xr.DataArray(d["X"].copy()), data2d)


def _assert_impute(astype, data):
    warnings.simplefilter("ignore", RuntimeWarning)
    Xorig = astype(data)
    X = copy(Xorig)
    assert_allclose(mean_impute(X), data["R"])
    assert_allclose(X, data["X"])
    assert isinstance(X, type(Xorig))

    Xorig = astype(data)
    X = copy(Xorig)
    assert_allclose(mean_impute(X, axis=0), data["Rt"])
    assert_allclose(X, data["X"])
    assert isinstance(X, type(Xorig))


def _assert_impute_inplace(astype, data):
    warnings.simplefilter("ignore", RuntimeWarning)
    Xorig = astype(data)
    X = copy(Xorig)
    assert_allclose(mean_impute(X, inplace=True), data["R"])
    assert_allclose(X, data["R"])
    assert isinstance(X, type(Xorig))

    Xorig = astype(data)
    X = copy(Xorig)
    assert_allclose(mean_impute(X, axis=0, inplace=True), data["Rt"])
    assert_allclose(X, data["Rt"])
    assert isinstance(X, type(Xorig))
