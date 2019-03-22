import sys
from limix._display import session_line

from .._data import conform_dataset
from .._data import asarray as _asarray
from .._display import session_block
from ._result import ScanResultFactory
from .._bits.xarray import is_dataarray


def st_iscan(G, y, E1, idx=None, K=None, M=None, E0=None, verbose=True):
    from glimix_core.lmm import LMM
    from numpy import ones, asarray
    from numpy_sugar.linalg import economic_qs
    from xarray import concat
    from numpy_sugar import is_all_finite

    if E0 is None:
        E0 = ones([G.shape[0], 1])

    E0 = _asarray(E0, "inter0", ["sample", "inter"])
    E1 = _asarray(E1, "inter1", ["sample", "inter"])
    E01 = concat([E0, E1], dim="inter")

    with session_block("qtl analysis", disable=not verbose):

        with session_line("Normalising input... ", disable=not verbose):
            data = conform_dataset(y, M, G=G, K=K)

        y = data["y"]
        M = data["M"]
        G = data["G"]
        K = data["K"]

        if idx is None:
            idx = range(G.shape[1])

        if not is_all_finite(data["y"]):
            raise ValueError("Outcome must have finite values only.")

        if not is_all_finite(data["M"]):
            raise ValueError("Covariates must have finite values only.")

        if data["K"] is not None:
            if not is_all_finite(data["K"]):
                raise ValueError("Covariate matrix must have finite values only.")
            QS = economic_qs(data["K"])
        else:
            QS = None

        lmm = LMM(y.values, M.values, QS)
        lmm.fit(verbose=verbose)
        sys.stdout.flush()

        r = ST_IScanResultFactory(M.covariate, G.candidate, E1.inter, E0.inter)
        r.set_null(lmm.lml(), lmm.beta, lmm.v1, lmm.v0)
        flmm = lmm.get_fast_scanner()
        for i in idx:

            i = _2d_sel(i)
            g = G[:, i]

            g0 = _cartesian_col(E0, g)
            g1 = _cartesian_col(E01, g)

            lml0, effs0 = flmm.scan(asarray(g0, float))
            lml1, effs1 = flmm.scan(asarray(g1, float))

            r.add_test(i, effs0, lml0, effs1, lml1)

        r = r.create()
        if verbose:
            print(r)

        return r


def _repeat(a, repeats, axis=None):
    try:
        import dask.array as da

        return da.repeat(a, repeats, axis)
    except (TypeError, AttributeError):
        import numpy as np

        return np.repeat(a, repeats, axis)


def _cartesian_col(a, b):
    from dask.array import tile

    na = a.shape[1]
    nb = b.shape[1]

    a = tile(a, nb)
    b = _repeat(b, na, axis=1)

    if is_dataarray(a):
        a = a.data

    if is_dataarray(b):
        b = b.data

    return a * b


def _2d_sel(idx):
    from collections.abc import Iterable

    if not isinstance(idx, (slice, Iterable)):
        return [idx]

    return idx
