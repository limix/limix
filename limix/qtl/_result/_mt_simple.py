from limix._bits import unvec
from limix._cache import cache


class MTSimpleModelResult:
    def __init__(self, lik, traits, covariates, lml, beta, beta_se, C0, C1):
        from numpy import asarray, atleast_1d, atleast_2d

        self._lik = lik
        self._traits = asarray(atleast_1d(traits), str)
        self._covariates = asarray(atleast_1d(covariates), str)
        self._lml = float(lml)
        self._beta = atleast_1d(asarray(beta, float).T).T
        self._beta_se = atleast_1d(asarray(beta_se, float).T).T
        self._C0 = atleast_2d(asarray(C0, float))
        self._C1 = atleast_2d(asarray(C1, float))

    @property
    def traits(self):
        return self._traits

    @property
    def likelihood(self):
        return self._lik

    @property
    def lml(self):
        return self._lml

    @property
    def effsizes(self):
        return self._dataframes["effsizes"]

    @property
    def variances(self):
        return self._dataframes["variances"]

    @property
    @cache
    def _dataframes(self):
        from pandas import DataFrame

        effsizes = []
        B = unvec(self._beta, (len(self._covariates), -1))
        Bvar = unvec(self._beta_se, (len(self._covariates), -1))
        for i, trait in enumerate(self._traits):
            for j, c in enumerate(self._covariates):
                effsizes.append([trait, c, B[j, i], Bvar[j, i]])

        columns = ["trait", "covariate", "effsize", "effsize_se"]
        df0 = DataFrame(effsizes, columns=columns)

        variances = []
        for i, trait0 in enumerate(self._traits):
            for j, trait1 in enumerate(self._traits):
                variances.append([trait0, trait1, self._C0[i, j], self._C1[i, j]])

        columns = ["trait0", "trait1", "fore_covariance", "back_covariance"]
        df1 = DataFrame(variances, columns=columns)

        return {"effsizes": df0, "variances": df1}
