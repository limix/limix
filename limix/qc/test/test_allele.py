import dask.array as da
import xarray as xr
from numpy import array_equal, ndarray
from numpy.random import RandomState
from numpy.testing import assert_, assert_allclose
from pandas import DataFrame, Series

from limix.qc import compute_maf


def test_compute_maf_numpy():
    random = RandomState(0)
    X = random.randint(0, 3, size=(100, 10))

    maf = compute_maf(X)
    assert_allclose(maf, [0.49, 0.49, 0.445, 0.495, 0.5, 0.45, 0.48, 0.48, 0.47, 0.435])


def test_compute_maf_dataframe():
    random = RandomState(0)

    X = random.randint(0, 3, size=(100, 10))
    columns = [f"snp{i}" for i in range(X.shape[1])]
    maf = compute_maf(DataFrame(X, columns=columns))

    assert_(isinstance(maf, Series))
    assert_(maf.name == "maf")
    assert_(array_equal(maf.index, columns))
    assert_allclose(maf, [0.49, 0.49, 0.445, 0.495, 0.5, 0.45, 0.48, 0.48, 0.47, 0.435])


def test_compute_maf_dask_array():
    random = RandomState(0)

    X = da.from_array(random.randint(0, 3, size=(100, 10)), chunks=2)
    maf = compute_maf(X)

    assert_(isinstance(maf, ndarray))
    assert_allclose(maf, [0.49, 0.49, 0.445, 0.495, 0.5, 0.45, 0.48, 0.48, 0.47, 0.435])


def test_compute_maf_dataarray():
    random = RandomState(0)

    X = random.randint(0, 3, size=(100, 10))
    samples = [f"snp{i}" for i in range(X.shape[0])]
    candidates = [f"snp{i}" for i in range(X.shape[1])]
    X = xr.DataArray(
        X,
        dims=["sample", "candidate"],
        coords={"sample": samples, "candidate": candidates},
    )
    maf = compute_maf(X)

    assert_(isinstance(maf, xr.DataArray))
    assert_(maf.name == "maf")
    assert_(array_equal(maf.candidate, candidates))
    assert_allclose(maf, [0.49, 0.49, 0.445, 0.495, 0.5, 0.45, 0.48, 0.48, 0.47, 0.435])
