import pytest
from numpy import asarray, nan
from numpy.random import RandomState

from limix.qc import mean_impute

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

    samples = [f"sample{i}" for i in range(X.shape[0])]
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
    samples = [f"sample{i}" for i in range(X.shape[0])]
    return {"X": X, "R": R, "Rt": Rt, "samples": samples, "candidates": None}


def test_impute_ndarray_1d(data1d):
    assert_ndarray_1d(mean_impute, data1d)


def test_impute_ndarray_2d(data2d):
    assert_ndarray_2d(mean_impute, data2d)


def test_impute_pandas_series(data1d):
    assert_pandas_series(mean_impute, data1d)


def test_impute_pandas_dataframe(data2d):
    assert_pandas_dataframe(mean_impute, data2d)


def test_impute_dask_array(data2d):
    assert_dask_array(mean_impute, data2d)


def test_impute_dask_dataframe(data2d):
    assert_dask_dataframe(mean_impute, data2d)


def test_impute_xarray_dataarray(data2d):
    assert_xarray_dataarray(mean_impute, data2d)
