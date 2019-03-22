import sys

from limix._display import session_line

from .._data import asarray as _asarray, conform_dataset
from .._display import session_block
from ._result import ScanResultFactory


def mt_scan(G, Y, idx=None, K=None, M=None, A=None, A1=None, A0=None, verbose=True):
    from glimix_core.lmm import RKron2Sum
    from numpy_sugar.linalg import economic_qs, ddot
    from xarray import concat
    from numpy_sugar import is_all_finite
    from numpy import eye, sqrt, asarray, empty, concatenate

    with session_block("qtl analysis", disable=not verbose):

        with session_line("Normalising input... ", disable=not verbose):

            data = conform_dataset(Y, M, G=G, K=K)

        Y = data["y"]
        M = data["M"]
        G = data["G"]
        K = data["K"]

        if A is None:
            A = eye(Y.shape[1])

        if A1 is None:
            A1 = eye(Y.shape[1])

        if A0 is None:
            A0 = empty((Y.shape[1], 0))

        A0 = _asarray(A0, "env0", ["sample", "env"])
        A1 = _asarray(A1, "env1", ["sample", "env"])
        A01 = concat([A0, A1], dim="env")

        if idx is None:
            idx = range(G.shape[1])

        if not is_all_finite(data["y"]):
            raise ValueError("Outcome must have finite values only.")

        if not is_all_finite(data["M"]):
            raise ValueError("Covariates must have finite values only.")

        if K is not None:
            if not is_all_finite(K):
                raise ValueError("Covariate matrix must have finite values only.")
            QS = economic_qs(K)
        else:
            QS = None

        X = ddot(QS[0][0], sqrt(QS[1]))
        lmm = RKron2Sum(Y.values, A, M.values, X)
        lmm.fit(verbose=verbose)
        sys.stdout.flush()

        flmm = lmm.get_fast_scanner()
        r = MT_ScanResultFactory(Y.trait, M.covariate, G.candidate, A1.env, A0.env)
        r.set_null(flmm.null_lml(), flmm.null_effsizes())
        for i in idx:

            i = _2d_sel(i)
            g = asarray(G[:, i], float)

            lml0, _, effs0_g = flmm.scan(A0, g)
            lml1, _, effs1_g = flmm.scan(A01, g)

            r.add_test(i, effs0_g, lml0, effs1_g, lml1)

        r = r.create()
        if verbose:
            print(r)

        return r


def _2d_sel(idx):
    from collections.abc import Iterable

    if not isinstance(idx, (slice, Iterable)):
        return [idx]

    return idx
