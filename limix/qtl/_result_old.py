from functools import lru_cache

from limix.stats import effsizes_se, lrt_pvalues


class ST_ScanResultFactory:
    def __init__(self, lik, covariates, candidates):
        from numpy import asarray

        self._lik = lik
        self._covariates = asarray(covariates, str)
        self._candidates = asarray(candidates, str)
        self._tests = []
        self._null_lml = None
        self._effsizes = None
        self._back_var = None
        self._fore_var = 0.0

    def set_null(self, null_lml, effsizes, back_var, fore_var=0.0):
        self._null_lml = null_lml
        self._effsizes = effsizes
        self._back_var = back_var
        self._fore_var = fore_var

    def add_test(self, cand_idx, alt_effsizes, alt_lml):
        from numpy import atleast_1d

        if not isinstance(cand_idx, slice):
            cand_idx = atleast_1d(cand_idx).ravel()

        alt_effsizes = _make_sure_iterabble(alt_effsizes)

        alt = {"lml": alt_lml, "effsizes": alt_effsizes}

        self._tests.append({"idx": cand_idx, "alt": alt})

    def create(self):
        return ST_ScanResult(
            self._tests,
            self._lik,
            self._covariates,
            self._candidates,
            self._null_lml,
            self._effsizes,
            {"background": self._back_var, "foreground": self._fore_var},
        )


class ST_ScanResult:
    def __init__(
        self, tests, lik, covariates, candidates, null_lml, effsizes, variances
    ):
        from numpy import asarray

        self._tests = tests
        self._lik = lik
        self._covariates = asarray(covariates, str)
        self._candidates = asarray(candidates, str)
        self._null_lml = null_lml
        self._effsizes = effsizes
        self._back_var = variances["background"]
        self._fore_var = variances["foreground"]

    @property
    def foreground_variance(self):
        return self._fore_var

    @property
    def background_variance(self):
        return self._back_var

    @property
    def stats(self):
        return self._dataframes["stats"].set_index("test")

    @property
    def alt_effsizes(self):
        return self._dataframes["effsizes"]["alt"]

    @property
    def covariate_effsizes(self):
        from pandas import DataFrame

        data = zip(self._covariates, self._effsizes)
        return DataFrame(data, columns=["covariate", "effsizes"]).set_index("covariate")

    @property
    @lru_cache(maxsize=None)
    def _dataframes(self):
        from pandas import DataFrame

        stats = []
        alt = []
        for i, test in enumerate(self._tests):
            idx = test["idx"]
            candidates = self._candidates[idx]

            alt_lml = test["alt"]["lml"]
            effsizes = test["alt"]["effsizes"]
            dof = len(effsizes)

            stats.append([i, self._null_lml, alt_lml, dof])
            alt += [[i, c, e] for c, e in zip(candidates, list(effsizes))]

        columns = ["test", "null lml", "alt lml", "dof"]
        stats = DataFrame(stats, columns=columns)

        columns = ["test", "candidate", "effsize"]
        alt = DataFrame(alt, columns=columns)

        stats["pvalue"] = lrt_pvalues(stats["null lml"], stats["alt lml"], stats["dof"])

        test = alt["test"]
        df = stats.set_index("test")
        alt["effsize se"] = effsizes_se(
            alt["effsize"], df.iloc[test]["pvalue"], df.iloc[test]["dof"]
        )

        stats = stats[["test", "null lml", "alt lml", "pvalue", "dof"]]
        return {"stats": stats, "effsizes": {"alt": alt}}

    def __repr__(self):
        from textwrap import TextWrapper
        from numpy import percentile, asarray

        msg = "Null model\n"
        msg += "----------\n\n"
        v0 = self.background_variance
        v1 = self.foreground_variance

        if self._lik[0] == "normal":
            msg += f"  𝐲 ~ 𝓝(M𝜶, {v0:.2f}*K + {v1:.2f}*I)\n"
        else:
            msg += f"  𝐳 ~ 𝓝(M𝜶, {v0:.2f}*K + {v1:.2f}*I)\n"

        msg += _lik_formulae(self._lik[0])

        prefix = "  M = "
        s = " " * len(prefix)
        wrapper = TextWrapper(initial_indent=prefix, width=88, subsequent_indent=s)
        msg += wrapper.fill(str(self._covariates)) + "\n"

        prefix = "  𝜶 = "
        s = " " * len(prefix)
        wrapper = TextWrapper(initial_indent=prefix, width=88, subsequent_indent=s)

        effsizes = asarray(self.covariate_effsizes, float).ravel()
        stats = self.stats

        msg += wrapper.fill(str(effsizes)) + "\n"
        msg += "  Log marg. lik.: {}\n".format(self._null_lml)
        msg += "  Number of models: 1\n\n"
        msg += "Alt model\n"
        msg += "---------\n\n"

        if self._lik[0] == "normal":
            msg += f"  𝐲 ~ 𝓝(M𝜶 + Gᵢ, {v0:.2f}*K + {v1:.2f}*I)\n"
        else:
            msg += f"  𝐳 ~ 𝓝(M𝜶 + Gᵢ, {v0:.2f}*K + {v1:.2f}*I)\n"

        msg += _lik_formulae(self._lik[0])

        msg += "  Min. p-value: {}\n".format(min(stats["pvalue"]))
        msg += "  First perc. p-value: {}\n".format(percentile(stats["pvalue"], 1))
        msg += "  Max. log marg. lik.: {}\n".format(max(stats["alt lml"]))
        msg += "  99th perc. log marg. lik.: {}\n".format(
            percentile(stats["alt lml"], 99)
        )
        msg += "  Number of models: {}".format(stats.shape[0])

        return msg


class MT_ScanResultFactory:
    def __init__(self, traits, covariates, candidates, envs0, envs1):
        from numpy import asarray

        self._traits = asarray(traits, str)
        self._covariates = asarray(covariates, str)
        self._candidates = asarray(candidates, str)
        self._envs0 = asarray(envs0)
        self._envs1 = asarray(envs1)
        self._tests = []
        self._null_lml = None
        self._effsizes = None

    def set_null(self, null_lml, effsizes):
        self._null_lml = null_lml
        self._effsizes = effsizes

    def add_test(self, cand_idx, null_effsizes, null_lml, alt_effsizes, alt_lml):
        from numpy import atleast_1d

        if not isinstance(cand_idx, slice):
            cand_idx = atleast_1d(cand_idx).ravel()

        null_effsizes = _make_sure_iterabble(null_effsizes)
        alt_effsizes = _make_sure_iterabble(alt_effsizes)

        null = {"lml": null_lml, "effsizes": null_effsizes}
        alt = {"lml": alt_lml, "effsizes": alt_effsizes}

        self._tests.append({"idx": cand_idx, "null": null, "alt": alt})

    def create(self):
        return MT_ScanResult(
            self._traits,
            self._covariates,
            self._candidates,
            self._envs1,
            self._envs0,
            self._tests,
            self._null_lml,
            self._effsizes,
        )


class MT_ScanResult:
    def __init__(
        self, traits, covariates, candidates, envs0, envs1, tests, null_lml, effsizes
    ):
        from numpy import asarray

        self._traits = traits
        self._covariates = asarray(covariates, str)
        self._candidates = asarray(candidates, str)
        self._envs0 = asarray(envs0)
        self._envs1 = asarray(envs1)
        self._tests = tests
        self._null_lml = null_lml
        self._effsizes = effsizes

    @property
    def stats(self):
        return self._dataframes["stats"].set_index("test")

    @property
    def null_effsizes(self):
        return self._dataframes["effsizes"]["null"]

    @property
    def alt_effsizes(self):
        return self._dataframes["effsizes"]["alt"]

    @property
    def covariate_effsizes(self):
        from pandas import DataFrame

        data = zip(self._covariates, self._effsizes)
        return DataFrame(data, columns=["covariate", "effsizes"]).set_index("covariate")

    @property
    @lru_cache(maxsize=None)
    def _dataframes(self):
        from pandas import DataFrame
        from itertools import product

        envs0 = list(self._envs0)
        envs1 = list(self._envs1)

        stats = []
        null = []
        alt = []
        for i, test in enumerate(self._tests):
            idx = test["idx"]
            candidates = list(self._candidates[idx])

            effsizes = list(test["null"]["effsizes"])

            for ci, eff in zip(list(product(candidates, envs0)), effsizes):
                null.append([i, ci[0], "env0_" + str(ci[1]), eff])

            effsizes = list(test["alt"]["effsizes"])
            for c in candidates:
                for j, env in enumerate(envs1):
                    alt.append([i, c, "env1_" + str(env), effsizes[j]])
                effsizes = effsizes[len(envs1) :]

            dof = len(test["alt"]["effsizes"]) - len(test["null"]["effsizes"])
            null_dof = len(test["null"]["effsizes"])

            stats.append([i, test["null"]["lml"], test["alt"]["lml"], dof, null_dof])

        columns = ["test", "null lml", "alt lml", "dof", "null dof"]
        stats = DataFrame(stats, columns=columns)

        columns = ["test", "candidate", "env", "effsize"]
        null = DataFrame(null, columns=columns)

        columns = ["test", "candidate", "env", "effsize"]
        alt = DataFrame(alt, columns=columns)

        stats["pvalue"] = lrt_pvalues(stats["null lml"], stats["alt lml"], stats["dof"])
        stats["null pvalue"] = lrt_pvalues(
            self._null_lml, stats["null lml"], stats["null dof"]
        )

        test = null["test"]
        df = stats.set_index("test")
        null["effsize se"] = effsizes_se(
            null["effsize"], df.iloc[test]["null pvalue"], df.iloc[test]["null dof"]
        )
        del stats["null pvalue"]
        del stats["null dof"]

        test = alt["test"]
        df = stats.set_index("test")
        alt["effsize se"] = effsizes_se(
            alt["effsize"], df.iloc[test]["pvalue"], df.iloc[test]["dof"]
        )

        return {"stats": stats, "effsizes": {"null": null, "alt": alt}}

    def __repr__(self):
        from textwrap import TextWrapper
        from numpy import percentile, asarray

        msg = "Null model\n"
        msg += "----------\n\n"
        v0 = self.background_variance
        v1 = self.foreground_variance
        msg += f"  y ~ 𝓝(M𝜶 + (E₀⊙Gᵢ)𝞫, {v0:.2f}*K + {v1:.2f}*I)\n"

        prefix = "  M = "
        s = " " * len(prefix)
        wrapper = TextWrapper(initial_indent=prefix, width=88, subsequent_indent=s)
        msg += wrapper.fill(str(self._covariates)) + "\n"

        prefix = "  𝜶 = "
        s = " " * len(prefix)
        wrapper = TextWrapper(initial_indent=prefix, width=88, subsequent_indent=s)

        effsizes = asarray(self.covariate_effsizes, float).ravel()
        stats = self.stats

        msg += wrapper.fill(str(effsizes)) + "\n"
        msg += "  Log marg. lik.: {}\n".format(self._null_lml)
        msg += "  Number of models: {}\n\n".format(stats.shape[0])
        msg += "Alt model\n"
        msg += "---------\n\n"
        msg += f"  y ~ 𝓝(M𝜶 + (E₀⊙Gᵢ)𝞫₀ + (E₁⊙Gᵢ)𝞫₁, {v0:.2f}*K + {v1:.2f}*I)\n"

        msg += "  Min. p-value: {}\n".format(min(stats["pvalue"]))
        msg += "  First perc. p-value: {}\n".format(percentile(stats["pvalue"], 1))
        msg += "  Max. log marg. lik.: {}\n".format(max(stats["alt lml"]))
        msg += "  99th perc. log marg. lik.: {}\n".format(
            percentile(stats["alt lml"], 99)
        )
        msg += "  Number of models: {}".format(stats.shape[0])

        return msg


class ST_IScanResultFactory:
    def __init__(self, covariates, candidates, inters1, inters0):
        from numpy import asarray

        self._covariates = asarray(covariates, str)
        self._candidates = asarray(candidates, str)
        self._inters0 = asarray(inters0)
        self._inters1 = asarray(inters1)
        self._tests = []
        self._null_lml = None
        self._effsizes = None
        self._back_var = None
        self._fore_var = 0.0

    def set_null(self, null_lml, effsizes, back_var, fore_var=0.0):
        self._null_lml = null_lml
        self._effsizes = effsizes
        self._back_var = back_var
        self._fore_var = fore_var

    def add_test(self, cand_idx, null_effsizes, null_lml, alt_effsizes, alt_lml):
        from numpy import atleast_1d

        if not isinstance(cand_idx, slice):
            cand_idx = atleast_1d(cand_idx).ravel()

        null_effsizes = _make_sure_iterabble(null_effsizes)
        alt_effsizes = _make_sure_iterabble(alt_effsizes)

        null = {"lml": null_lml, "effsizes": null_effsizes}
        alt = {"lml": alt_lml, "effsizes": alt_effsizes}

        self._tests.append({"idx": cand_idx, "null": null, "alt": alt})

    def create(self):
        return ST_IScanResult(
            self._covariates,
            self._candidates,
            self._inters1,
            self._inters0,
            self._tests,
            self._null_lml,
            self._effsizes,
            {"background": self._back_var, "foreground": self._fore_var},
        )


class ST_IScanResult:
    def __init__(
        self,
        covariates,
        candidates,
        inters1,
        inters0,
        tests,
        null_lml,
        effsizes,
        variances,
    ):
        self._covariates = covariates
        self._candidates = candidates
        self._inters0 = inters0
        self._inters1 = inters1
        self._tests = tests
        self._null_lml = null_lml
        self._effsizes = effsizes
        self._back_var = variances["background"]
        self._fore_var = variances["foreground"]

    @property
    def foreground_variance(self):
        return self._fore_var

    @property
    def background_variance(self):
        return self._back_var

    @property
    def stats(self):
        return self._dataframes["stats"].set_index("test")

    @property
    def null_effsizes(self):
        return self._dataframes["effsizes"]["null"]

    @property
    def alt_effsizes(self):
        return self._dataframes["effsizes"]["alt"]

    @property
    def covariate_effsizes(self):
        from pandas import DataFrame

        data = zip(self._covariates, self._effsizes)
        return DataFrame(data, columns=["covariate", "effsizes"]).set_index("covariate")

    @property
    @lru_cache(maxsize=None)
    def _dataframes(self):
        from pandas import DataFrame
        from itertools import product

        inters0 = list(self._inters0)
        inters1 = list(self._inters1)

        stats = []
        null = []
        alt = []
        for i, test in enumerate(self._tests):
            idx = test["idx"]
            candidates = list(self._candidates[idx])

            effsizes = list(test["null"]["effsizes"])

            for ci, eff in zip(list(product(candidates, inters0)), effsizes):
                null.append([i, ci[0], "inter0_" + str(ci[1]), eff])

            effsizes = list(test["alt"]["effsizes"])
            for c in candidates:

                for j, inter in enumerate(inters0):
                    alt.append([i, c, "inter0_" + str(inter), effsizes[j]])
                effsizes = effsizes[len(inters0) :]

                for j, inter in enumerate(inters1):
                    alt.append([i, c, "inter1_" + str(inter), effsizes[j]])
                effsizes = effsizes[len(inters1) :]

            dof = len(test["alt"]["effsizes"]) - len(test["null"]["effsizes"])
            null_dof = len(test["null"]["effsizes"])

            stats.append([i, test["null"]["lml"], test["alt"]["lml"], dof, null_dof])

        columns = ["test", "null lml", "alt lml", "dof", "null dof"]
        stats = DataFrame(stats, columns=columns)

        columns = ["test", "candidate", "inter", "effsize"]
        null = DataFrame(null, columns=columns)

        columns = ["test", "candidate", "inter", "effsize"]
        alt = DataFrame(alt, columns=columns)

        stats["pvalue"] = lrt_pvalues(stats["null lml"], stats["alt lml"], stats["dof"])
        stats["null pvalue"] = lrt_pvalues(
            self._null_lml, stats["null lml"], stats["null dof"]
        )

        test = null["test"]
        df = stats.set_index("test")
        null["effsize se"] = effsizes_se(
            null["effsize"], df.iloc[test]["null pvalue"], df.iloc[test]["null dof"]
        )
        del stats["null pvalue"]
        del stats["null dof"]

        test = alt["test"]
        df = stats.set_index("test")
        alt["effsize se"] = effsizes_se(
            alt["effsize"], df.iloc[test]["pvalue"], df.iloc[test]["dof"]
        )

        return {"stats": stats, "effsizes": {"null": null, "alt": alt}}

    def __repr__(self):
        from textwrap import TextWrapper
        from numpy import percentile, asarray

        msg = "Null model\n"
        msg += "----------\n\n"
        v0 = self.background_variance
        v1 = self.foreground_variance
        msg += f"  y ~ 𝓝(M𝜶 + (E₀⊙Gᵢ)𝞫, {v0:.2f}*K + {v1:.2f}*I)\n"

        prefix = "  M = "
        s = " " * len(prefix)
        wrapper = TextWrapper(initial_indent=prefix, width=88, subsequent_indent=s)
        msg += wrapper.fill(str(self._covariates)) + "\n"

        prefix = "  𝜶 = "
        s = " " * len(prefix)
        wrapper = TextWrapper(initial_indent=prefix, width=88, subsequent_indent=s)

        effsizes = asarray(self.covariate_effsizes, float).ravel()
        stats = self.stats

        msg += wrapper.fill(str(effsizes)) + "\n"
        msg += "  Log marg. lik.: {}\n".format(self._null_lml)
        msg += "  Number of models: {}\n\n".format(stats.shape[0])
        msg += "Alt model\n"
        msg += "---------\n\n"
        msg += f"  y ~ 𝓝(M𝜶 + (E₀⊙Gᵢ)𝞫₀ + (E₁⊙Gᵢ)𝞫₁, {v0:.2f}*K + {v1:.2f}*I)\n"

        msg += "  Min. p-value: {}\n".format(min(stats["pvalue"]))
        msg += "  First perc. p-value: {}\n".format(percentile(stats["pvalue"], 1))
        msg += "  Max. log marg. lik.: {}\n".format(max(stats["alt lml"]))
        msg += "  99th perc. log marg. lik.: {}\n".format(
            percentile(stats["alt lml"], 99)
        )
        msg += "  Number of models: {}".format(stats.shape[0])

        return msg


def _make_sure_iterabble(x):
    from collections.abc import Iterable

    if not isinstance(x, Iterable):
        return [x]

    return x


def _lik_formulae(lik):
    msg = ""

    if lik == "bernoulli":
        msg += f"  yᵢ ~ Bern(μᵢ=g(zᵢ)), where g(x)=1/(1+e⁻ˣ)\n"
    elif lik == "probit":
        msg += f"  yᵢ ~ Bern(μᵢ=g(zᵢ)), where g(x)=Φ(x)\n"
    elif lik == "binomial":
        msg += f"  yᵢ ~ Binom(μᵢ=g(zᵢ), nᵢ), where g(x)=1/(1+e⁻ˣ)\n"
    elif lik == "poisson":
        msg += f"  yᵢ ~ Poisson(λᵢ=g(zᵢ)), where g(x)=eˣ\n"

    return msg
