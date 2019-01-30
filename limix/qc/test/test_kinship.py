import dask.array as da
import xarray as xr
from numpy import ndarray, zeros
from numpy.random import RandomState
from numpy.testing import assert_, assert_allclose
from pandas import DataFrame

from limix.qc import normalise_covariance


def test_qc_kinship_numpy():
    random = RandomState(0)
    X = random.randn(3, 5)
    K = X.dot(X.T)
    K1 = zeros((3, 3))

    K0 = normalise_covariance(K)
    K2 = normalise_covariance(K, out=K1)

    Kf = [
        [2.5990890007787586, -0.1951278087849671, 0.5472860002747189],
        [-0.1951278087849671, 0.4202620710126438, 0.2642930556468809],
        [0.5472860002747189, 0.2642930556468809, 0.5971001753452302],
    ]

    assert_allclose(K0, Kf)
    assert_allclose(K0, K1)
    assert_allclose(K0, K2)
    assert_(K2 is K1)


def test_qc_kinship_dataframe():
    random = RandomState(0)
    X = random.randn(3, 5)
    K = X.dot(X.T)
    K = DataFrame(K)

    K1 = DataFrame(zeros((3, 3)))

    K0 = normalise_covariance(K)
    K2 = normalise_covariance(K, out=K1)

    Kf = [
        [2.5990890007787586, -0.1951278087849671, 0.5472860002747189],
        [-0.1951278087849671, 0.4202620710126438, 0.2642930556468809],
        [0.5472860002747189, 0.2642930556468809, 0.5971001753452302],
    ]

    assert_allclose(K0, Kf)
    assert_(isinstance(K0, DataFrame))

    assert_allclose(K0, K1)
    assert_(isinstance(K1, DataFrame))

    assert_allclose(K0, K2)
    assert_(K2 is K1)


def test_qc_kinship_dask_array():
    random = RandomState(0)
    X = random.randn(3, 5)
    K = X.dot(X.T)
    K = da.from_array(K, chunks=2)

    K1 = zeros((3, 3))

    K0 = normalise_covariance(K)
    K2 = normalise_covariance(K, out=K1)

    Kf = [
        [2.5990890007787586, -0.1951278087849671, 0.5472860002747189],
        [-0.1951278087849671, 0.4202620710126438, 0.2642930556468809],
        [0.5472860002747189, 0.2642930556468809, 0.5971001753452302],
    ]

    assert_allclose(K0, Kf)
    assert_(isinstance(K0, da.Array))

    assert_allclose(K0, K1)
    assert_(isinstance(K1, ndarray))
    assert_(isinstance(K2, ndarray))

    assert_allclose(K0, K2)
    assert_(K2 is K1)


def test_qc_kinship_dataarray():
    random = RandomState(0)
    X = random.randn(3, 5)
    K = X.dot(X.T)
    K = da.from_array(K, chunks=2)
    K = xr.DataArray(K)

    K1 = zeros((3, 3))

    K0 = normalise_covariance(K)
    K2 = normalise_covariance(K, out=K1)

    Kf = [
        [2.5990890007787586, -0.1951278087849671, 0.5472860002747189],
        [-0.1951278087849671, 0.4202620710126438, 0.2642930556468809],
        [0.5472860002747189, 0.2642930556468809, 0.5971001753452302],
    ]

    assert_allclose(K0, Kf)
    assert_(isinstance(K0, xr.DataArray))

    assert_allclose(K0, K1)
    assert_(isinstance(K1, ndarray))
    assert_(isinstance(K2, ndarray))

    assert_allclose(K0, K2)
    assert_(K2 is K1)
