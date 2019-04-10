from ._result import ScanResult
from ._simple import SModelResult


class STScanResultFactory:
    def __init__(self, lik, trait, covariates, candidates, lml, beta, beta_se, v0, v1):
        from numpy import asarray, atleast_1d

        self._h0 = SModelResult(lik, trait, covariates, lml, beta, beta_se, v0, v1)
        self._tests = []
        self._trait = str(trait)
        self._covariates = asarray(atleast_1d(covariates), str)
        self._candidates = asarray(atleast_1d(candidates), str)

    def add_test(self, cand_idx, h1):
        from numpy import atleast_1d, atleast_2d, asarray

        if not isinstance(cand_idx, slice):
            cand_idx = asarray(atleast_1d(cand_idx).ravel(), int)

        def _1d_shape(x):
            x = asarray(x, float)
            x = atleast_1d(x.T).T
            return x

        def _normalize(h):
            return {
                "lml": float(h["lml"]),
                "covariate_effsizes": _1d_shape(h["covariate_effsizes"]),
                "candidate_effsizes": _1d_shape(h["candidate_effsizes"]),
                "covariate_effsizes_se": _1d_shape(h["covariate_effsizes_se"]),
                "candidate_effsizes_se": _1d_shape(h["candidate_effsizes_se"]),
                "scale": float(h["scale"]),
            }

        self._tests.append({"idx": cand_idx, "h1": _normalize(h1)})

    def create(self):
        return ScanResult(
            self._tests, self._trait, self._covariates, self._candidates, self._h0
        )
