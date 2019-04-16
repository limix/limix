import dask.array as da
import dask.dataframe as dd
import pandas as pd
import xarray as xr
from numpy.testing import assert_, assert_allclose, assert_array_equal


def assert_ndarray_1d(func, data1d):
    _assert_mat_proc(func, data1d, lambda X, *_: X.copy())
    _assert_mat_proc_inplace(func, data1d, lambda X, *_: X.copy())


def assert_ndarray_2d(func, data2d):
    _assert_mat_proc(func, data2d, lambda X, *_: X.copy())
    _assert_mat_proc_inplace(func, data2d, lambda X, *_: X.copy())


def assert_pandas_series(func, data1d):
    def get(X, samples, *_):
        return pd.Series(X, index=samples)

    _assert_mat_proc(func, data1d, get)
    _assert_mat_proc_inplace(func, data1d, get)


def assert_pandas_dataframe(func, data2d):
    def get(X, samples, candidates):
        return pd.DataFrame(X, index=samples, columns=candidates)

    _assert_mat_proc(func, data2d, get)
    _assert_mat_proc_inplace(func, data2d, get)


def assert_dask_array(func, data2d, test_inplace=True):
    def get(X, *_):
        return da.from_array(X.copy(), chunks=2)

    _assert_mat_proc(func, data2d, get)
    if test_inplace:
        _assert_mat_proc_inplace(func, data2d, get)


def assert_dask_dataframe(func, data2d):
    def get(X, samples, candidates):
        df = pd.DataFrame(X.copy(), index=samples, columns=candidates)
        df = dd.from_pandas(df, npartitions=2)
        return df

    _assert_mat_proc(func, data2d, get)


def assert_xarray_dataarray(func, data2d, test_inplace=True):
    def get(X, samples, candidates):
        X = xr.DataArray(
            X.copy(),
            coords={"sample": samples, "candidate": candidates},
            dims=["sample", "candidate"],
        )
        return X

    _assert_mat_proc(func, data2d, get)
    if test_inplace:
        _assert_mat_proc_inplace(func, data2d, get)


def _assert_mat_proc(func, data, astype):
    def get():
        return astype(data["X"], data["samples"], data["candidates"]).copy()

    X = get()
    assert_allclose(func(X), data["R"])
    assert_allclose(X, get())
    assert_(isinstance(func(X), type(get())))
    _assert_attr_values(X, get())

    X = get()
    assert_allclose(func(X, axis=0), data["Rt"])
    assert_allclose(X, get())
    assert_(isinstance(func(X, axis=0), type(get())))
    _assert_attr_values(X, get())


def _assert_mat_proc_inplace(func, data, astype):
    def get():
        return astype(data["X"], data["samples"], data["candidates"]).copy()

    X = get().copy()
    assert_allclose(func(X, inplace=True), data["R"])
    assert_allclose(X, data["R"])
    assert_(isinstance(func(X, inplace=True), type(get())))
    _assert_attr_values(X, get())

    X = get().copy()
    assert_allclose(func(X, axis=0, inplace=True), data["Rt"])
    assert_allclose(X, data["Rt"])
    assert_(isinstance(func(X, axis=0, inplace=True), type(get())))
    _assert_attr_values(X, get())


def _assert_attr_values(X, R):
    for attr in ["index", "columns"]:
        if hasattr(X, attr):
            assert_array_equal(getattr(X, attr), getattr(R, attr))
