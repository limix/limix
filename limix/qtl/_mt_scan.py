import sys

from limix._display import session_line

from .._data import asarray as _asarray, conform_dataset
from .._display import session_block
from ._result import ScanResultFactory


def mt_scan(G, Y, idx=None, K=None, M=None, A=None, A0=None, A1=None, verbose=True):
    r"""
    H0: vec(Y) ~ N((AxM)vec(B), C0*K + C1*I)
    H1: vec(Y) ~ N((AxM)vec(beta_1i)+(A0xGi)vec(alpha_1i), s_1i(C0*K + C1*I))
    H2: vec(Y) ~ N((AxM)vec(beta_2i)+([A0 A1]xGi)vec(alpha_2i), s_2i(C0*K + C1*I))
    """
    from glimix_core.lmm import Kron2Sum
    from numpy_sugar.linalg import economic_qs, ddot
    from xarray import concat
    from ._assert import assert_finite
    from numpy import eye, sqrt, asarray, empty

    with session_block("qtl analysis", disable=not verbose):

        with session_line("Normalising input... ", disable=not verbose):

            data = conform_dataset(Y, M, G=G, K=K)

        Y = data["y"]
        M = data["M"]
        G = data["G"]
        K = data["K"]

        assert_finite(Y, M, K)

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

        if K is not None:
            QS = economic_qs(K)
        else:
            QS = None

        KG = ddot(QS[0][0], sqrt(QS[1]))
        lmm = Kron2Sum(Y.values, A, M.values, KG, restricted=False)
        lmm.fit(verbose=verbose)
        sys.stdout.flush()

        scanner = lmm.get_fast_scanner()
        r = ScanResultFactory(
            "normal",
            Y.trait,
            M.covariate,
            G.candidate,
            A0.env,
            A1.env,
            scanner.null_lml(),
            scanner.null_beta,
            sqrt(scanner.null_beta_covariance.diagonal()),
            lmm.C0,
            lmm.C1,
        )

        for i in idx:

            i = _2d_sel(i)
            g = asarray(G[:, i], float)

            r1 = scanner.scan(A0, g)
            r2 = scanner.scan(A01, g)

            h1 = _normalise_scan_names(r1)
            h2 = _normalise_scan_names(r2)
            r.add_test(i, h1, h2)

        r = r.create()
        if verbose:
            print(r)

        return r


def _normalise_scan_names(r):
    return {
        "lml": r["lml"],
        "covariate_effsizes": r["effsizes0"],
        "covariate_effsizes_se": r["effsizes0_se"],
        "candidate_effsizes": r["effsizes1"],
        "candidate_effsizes_se": r["effsizes1_se"],
        "scale": r["scale"],
    }


def _2d_sel(idx):
    from collections.abc import Iterable

    if not isinstance(idx, (slice, Iterable)):
        return [idx]

    return idx
