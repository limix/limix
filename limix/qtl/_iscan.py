import sys

from limix._display import session_line

from .._bits import unvec
from .._data import asarray as _asarray, assert_likelihood, conform_dataset
from .._display import session_block
from ._result import IScanResultFactory


def iscan(G, y, lik="normal", K=None, M=None, idx=None, E0=None, E1=None, verbose=True):
    from numpy_sugar.linalg import economic_qs
    from xarray import concat
    from ._assert import assert_finite
    from numpy import asarray, empty, ones

    if not isinstance(lik, (tuple, list)):
        lik = (lik,)

    lik_name = lik[0].lower()
    lik = (lik_name,) + lik[1:]
    assert_likelihood(lik_name)

    with session_block("QTL analysis", disable=not verbose):

        with session_line("Normalising input... ", disable=not verbose):

            data = conform_dataset(y, M, G=G, K=K)

        Y = data["y"]
        M = data["M"]
        G = data["G"]
        K = data["K"]

        assert_finite(y, M, K)
        nsamples = y.shape[0]

        if E1 is None:
            E1 = ones((nsamples, 1))

        if E0 is None:
            E0 = empty((nsamples, 0))

        E0 = _asarray(E0, "env0", ["sample", "env"])
        E1 = _asarray(E1, "env1", ["sample", "env"])
        E01 = concat([E0, E1], dim="env")

        if K is not None:
            QS = economic_qs(K)
        else:
            QS = None

        if lik_name == "normal":
            scanner, v0, v1 = _lmm(Y.values.ravel(), M.values, QS, verbose)
        else:
            scanner, v0, v1 = _glmm(Y.values.ravel(), lik, M.values, QS, verbose)

        r = IScanResultFactory(
            lik_name,
            Y.trait,
            M.covariate,
            G.candidate,
            E0.env,
            E1.env,
            scanner.null_lml,
            scanner.null_beta,
            scanner.null_beta_se,
            v0,
            v1,
        )

        if idx is None:

            if E0.shape[1] == 0:
                r1 = scanner.fast_scan(G, verbose)
            else:
                r1 = scanner.scan(G, E0, verbose)

            if E1.shape[1] == 0:
                r2 = r1
            else:
                r2 = scanner.scan(G, E01, verbose)

            for i in range(G.shape[1]):
                h1 = _normalise_scan_names({k: v[i] for k, v in r1.items()})
                h2 = _normalise_scan_names({k: v[i] for k, v in r2.items()})
                r.add_test(i, h1, h2)
        else:
            for i in idx:
                i = _2d_sel(i)
                g = asarray(G[:, i], float)

                r1 = scanner.scan(g, E0)
                r2 = scanner.scan(g, E01)

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


class ScannerWrapper:
    def __init__(self, scanner):
        self._scanner = scanner

    @property
    def null_lml(self):
        return self._scanner.null_lml()

    @property
    def null_beta(self):
        return self._scanner.null_beta

    @property
    def null_beta_se(self):
        from numpy import sqrt

        se = sqrt(self._scanner.null_beta_covariance.diagonal())
        return se

    def fast_scan(self, G, verbose):
        r = self._scanner.fast_scan(G, verbose=verbose)
        r["effsizes1"] = unvec(r["effsizes1"], (-1, G.shape[1])).T
        r["effsizes1_se"] = unvec(r["effsizes1_se"], (-1, G.shape[1])).T
        return r

    def scan(self, G, E):
        from .._bits import cdot

        r = self._scanner.scan(cdot(G, E))
        r["effsizes1"] = unvec(r["effsizes1"], (-1, G.shape[1])).T
        r["effsizes1_se"] = unvec(r["effsizes1_se"], (-1, G.shape[1])).T
        return r


def _lmm(y, M, QS, verbose):
    from glimix_core.lmm import LMM

    lmm = LMM(y, M, QS, restricted=False)
    lmm.fit(verbose=verbose)
    sys.stdout.flush()

    if QS is None:
        v0 = None
    else:
        v0 = lmm.v0
    v1 = lmm.v1
    scanner = ScannerWrapper(lmm.get_fast_scanner())

    return scanner, v0, v1


def _glmm(y, lik, M, QS, verbose):
    from glimix_core.glmm import GLMMExpFam, GLMMNormal

    glmm = GLMMExpFam(y.ravel(), lik, M, QS)

    glmm.fit(verbose=verbose)
    v0 = glmm.v0
    v1 = glmm.v1
    sys.stdout.flush()

    eta = glmm.site.eta
    tau = glmm.site.tau

    gnormal = GLMMNormal(eta, tau, M, QS)
    gnormal.fit(verbose=verbose)

    scanner = ScannerWrapper(gnormal.get_fast_scanner())

    return scanner, v0, v1
