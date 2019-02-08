import dask.array as da
import xarray as xr
from limix.qc import quantile_gaussianize
from numpy import array_equal, ndarray, asarray, nan
from numpy.random import RandomState
from numpy.testing import assert_, assert_allclose
from pandas import DataFrame, Series


def test_trans_gaussianize_numpy():
    random = RandomState(0)
    X = random.randint(0, 3, size=(100, 10))

    r0 = quantile_gaussianize(X[:, 0])
    assert_allclose(
        r0[:5],
        [
            -0.9221781782775866,
            -0.9221781782775866,
            -0.9221781782775866,
            0.024820650014443935,
            0.9608379310031608,
        ],
    )
    R = quantile_gaussianize(X)
    assert_allclose(R[:, 0], r0)
    assert_allclose(
        [R[0, -1], R[1, -1], R[2, -1], R[3, -1]],
        [
            -1.0644415649029644,
            -0.16202377853274785,
            -0.16202377853274785,
            0.8136568081151939,
        ],
    )
    X = X.astype(float)
    X[0, 0] = nan
    R = quantile_gaussianize(X)
    assert_allclose(R[:4, 0], [nan, -0.93458929107348, -0.93458929107348,
                               0.012533469508069276])

def test_trans_gaussianize_dataframe():
    random = RandomState(0)
    X = random.randint(0, 3, size=(100, 10))
    columns = [f"snp{i}" for i in range(X.shape[1])]
    X = DataFrame(X, columns=columns)

    R = quantile_gaussianize(X)
    X = R.to_numpy()
    assert_allclose(
        [X[0, -1], X[1, -1], X[2, -1], X[3, -1]],
        [
            -1.0644415649029644,
            -0.16202377853274785,
            -0.16202377853274785,
            0.8136568081151939,
        ],
    )
    assert_(isinstance(R, DataFrame))
    assert_(array_equal(R.columns, columns))



def test_compute_impute_dask_array():
    random = RandomState(0)
    X = random.randint(0, 3, size=(100, 10))
    columns = [f"snp{i}" for i in range(X.shape[1])]
    X = da.from_array(X, chunks=2)

    R = quantile_gaussianize(X)
    X = asarray(R, float)
    assert_allclose(
        [X[0, -1], X[1, -1], X[2, -1], X[3, -1]],
        [
            -1.0644415649029644,
            -0.16202377853274785,
            -0.16202377853274785,
            0.8136568081151939,
        ],
    )
    assert_(isinstance(R, da.Array))


# def test_compute_impute_dataarray():
#     random = RandomState(0)

#     X = random.randint(0, 3, size=(100, 10))
#     samples = [f"snp{i}" for i in range(X.shape[0])]
#     candidates = [f"snp{i}" for i in range(X.shape[1])]
#     X = xr.DataArray(
#         X,
#         dims=["sample", "candidate"],
#         coords={"sample": samples, "candidate": candidates},
#     )
#     maf = compute_maf(X)

#     assert_(isinstance(maf, xr.DataArray))
#     assert_(maf.name == "maf")
#     assert_(array_equal(maf.candidate, candidates))
#     assert_allclose(maf, [0.49, 0.49, 0.445, 0.495, 0.5, 0.45, 0.48, 0.48, 0.47, 0.435])
