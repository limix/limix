from ._mt_result import MTScanResult
from ._mt_simple import MTSimpleModelResult
from ._tuples import VariantResult, Result


class MTScanResultFactory:
    def __init__(
        self,
        lik,
        traits,
        covariates,
        candidates,
        envs0,
        envs1,
        lml,
        beta,
        beta_se,
        C0,
        C1,
    ):
        from numpy import asarray, atleast_1d

        self._h0 = MTSimpleModelResult(
            lik, traits, covariates, lml, beta, beta_se, C0, C1
        )
        self._tests = []
        self._traits = asarray(atleast_1d(traits), str)
        self._covariates = asarray(atleast_1d(covariates), str)
        self._candidates = asarray(atleast_1d(candidates), str)
        self._envs0 = asarray(atleast_1d(envs0), str)
        self._envs1 = asarray(atleast_1d(envs1), str)

    def add_test(self, cand_idx, h1, h2):
        from numpy import atleast_1d, atleast_2d, asarray

        if not isinstance(cand_idx, slice):
            cand_idx = asarray(atleast_1d(cand_idx).ravel(), int)

        def _2d_shape(x):
            x = asarray(x, float)
            x = atleast_2d(x.T).T
            return x

        def _normalize(h):
            return VariantResult(
                lml=float(h.lml),
                covariate_effsizes=_2d_shape(h.covariate_effsizes),
                candidate_effsizes=_2d_shape(h.candidate_effsizes),
                covariate_effsizes_se=_2d_shape(h.covariate_effsizes_se),
                candidate_effsizes_se=_2d_shape(h.candidate_effsizes_se),
                scale=float(h.scale),
            )

        if h1 is not None:
            h1 = _normalize(h1)
        self._tests.append(Result(idx=cand_idx, h1=h1, h2=_normalize(h2)))

    def create(self):
        return MTScanResult(
            self._tests,
            self._traits,
            self._covariates,
            self._candidates,
            self._h0,
            self._envs0,
            self._envs1,
        )
