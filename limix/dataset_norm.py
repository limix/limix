import xarray as xr
from pprint import pprint
from xarray import DataArray
from numpy import asarray
from numpy_sugar import is_all_finite
import collections


def _normalise_outcome(y):
    if hasattr(y, "ndim"):
        if y.ndim == 1:
            y = [y]
        elif y.ndim == 2:
            y = y.T
        else:
            raise ValueError("Unrecognized number of dimensions of the outcome array.")
    else:
        if isinstance(y, collections.Sequence):
            if not hasattr(y[0], "ndim") and not isinstance(y[0], collections.Sequence):
                y = [asarray(y, float)]
        else:
            y = [y]
    y = [DataArray(yi, encoding={"dtype": "float64"}) for yi in y]
    if not all(is_all_finite(yi) for yi in y):
        raise ValueError("There are non-finite values in the outcome.")

    return y


def normalise_dataset(y, lik, M=None, G=None, K=None):
    if K is None:
        K = []

    # if isinstance(y, (tuple, list)):
    #     y = asarray(y, float).T

    y = _normalise_outcome(y)
    return (y,)
    # if M is None:
    #     M = ones((y.shape[0], 1))
    #     if _isindexed(y):
    #         M = DataFrame(M, index=y.index)

    # arrs = [a for a in [y, G, M, K] if a is not None]

    # # Exit early if all the indices are the same.
    # if all([_isindexed(a) for a in arrs]):
    #     a0 = arrs[0]
    #     if all([len(_getindex(a)) == len(_getindex(a0)) for a in arrs[1:]]):
    #         if all([all(_getindex(a) == _getindex(a0)) for a in arrs[1:]]):
    #             return dict(y=y, M=M, G=G, K=K)

    # samples_index = _infer_samples_index(arrs)
    # arrs0 = [a for a in [G, M, K] if a is not None and _isindexed(a)]
    # ok = _check_duplicity(samples_index, arrs0)

    # if not ok:
    #     msg = "Duplicated indices are allowed only on the"
    #     msg += " outcome variable when the indices differ."
    #     raise ValueError(msg)

    # if len(samples_index) == 0:
    #     raise ValueError(
    #         "Could not infer an index for samples."
    #         " Please, check if the passed arrays are in order."
    #     )

    # y = _make_sure_phenotype_dataframe(y, samples_index)
    # samples_index = y.index

    # if G is not None:
    #     G = _make_sure_candidates_dataframe(G, samples_index)
    #     if not isinstance(G, (da.DataFrame, xr.DataArray)) and not npall(isfinite(G)):
    #         raise ValueError("Candidate values must be finite.")

    # M = _make_sure_covariates_dataframe(M, samples_index)

    # if K is not None:
    #     K = _make_sure_kinship_dataframe(K, samples_index)
    #     if not isinstance(K, da.DataFrame) and not npall(isfinite(K)):
    #         raise ValueError("Covariance values must be finite.")

    # y = normalise_extreme_values(y, lik)

    # return dict(y=y, M=M, G=G, K=K)


if __name__ == "__main__":
    import numpy as np
    import pandas as pd
    import xarray as xr

    sample_ids_0 = ["sample_{}".format(i) for i in range(5)]

    y = np.random.randn(5)
    pprint(normalise_dataset(y, "normal")[0])
    print()

    y = np.random.randn(5, 1)
    pprint(normalise_dataset(y, "normal")[0])
    print()

    y = np.random.randn(5, 2)
    pprint(normalise_dataset(y, "normal")[0])
    print()

    y = np.random.randn(5)
    y = pd.Series(y, index=sample_ids_0)
    pprint(normalise_dataset(y, "normal")[0])
    print()

    y = np.random.randn(5)
    y = pd.Series(y, index=sample_ids_0)
    pprint(normalise_dataset(y, "normal")[0])
    print()
