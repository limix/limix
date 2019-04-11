from limix._cache import cache


class STSimpleModelResult:
    def __init__(self, lik, trait, covariates, lml, beta, beta_se, v0, v1):
        from numpy import asarray, atleast_1d

        self._lik = lik
        self._trait = str(trait)
        self._covariates = asarray(atleast_1d(covariates), str)
        self._lml = float(lml)
        self._beta = asarray(beta, float).ravel()
        self._beta_se = asarray(beta_se, float).ravel()
        self._v0 = float(v0)
        self._v1 = float(v1)

    @property
    def trait(self):
        return self._trait

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
        for j, c in enumerate(self._covariates):
            effsizes.append([self._trait, c, self._beta[j], self._beta_se[j]])

        columns = ["trait", "covariate", "effsize", "effsize_se"]
        df0 = DataFrame(effsizes, columns=columns)

        variances = [[self._trait, self._v0, self._v1]]
        columns = ["trait", "fore_covariance", "back_covariance"]
        df1 = DataFrame(variances, columns=columns)

        return {"effsizes": df0, "variances": df1}
