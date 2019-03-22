from limix._cache import cache
from limix.stats import lrt_pvalues


class ScanResultFactory:
    def __init__(self, traits, covariates, candidates, envs0, envs1):
        from numpy import asarray, atleast_1d

        self._tests = []
        self._traits = asarray(atleast_1d(traits), str)
        self._covariates = asarray(atleast_1d(covariates), str)
        self._candidates = asarray(atleast_1d(candidates), str)
        self._envs0 = asarray(atleast_1d(envs0), str)
        self._envs1 = asarray(atleast_1d(envs1), str)
        self._h0 = {"lml": None, "effsizes": None, "C0": None, "C1": None}

    def set_h0(self, lml, effsizes, C0, C1):
        from numpy import asarray, atleast_2d

        self._h0["lml"] = lml
        self._h0["effsizes"] = asarray(atleast_2d(effsizes), float)
        self._h0["C0"] = asarray(atleast_2d(C0), float)
        self._h0["C1"] = asarray(atleast_2d(C1), float)

    def add_test(self, cand_idx, h1, h2):
        from numpy import atleast_1d, atleast_2d, asarray

        if not isinstance(cand_idx, slice):
            cand_idx = asarray(atleast_1d(cand_idx).ravel(), int)

        def _normalize(h):
            h = {
                "lml": float(h["lml"]),
                "covariate_effsizes": atleast_2d(h["covariate_effsizes"]),
                "candidate_effsizes": atleast_2d(h["candidate_effsizes"]),
            }
            h["covariate_effsizes"] = asarray(h["covariate_effsizes"], float)
            h["candidate_effsizes"] = asarray(h["candidate_effsizes"], float)
            return h

        h1 = _normalize(h1)
        h2 = _normalize(h2)

        self._tests.append({"idx": cand_idx, "h1": h1, "h2": h2})

    def create(self):
        return ScanResult(
            self._tests,
            self._traits,
            self._covariates,
            self._candidates,
            self._h0,
            self._envs0,
            self._envs1,
        )


class ScanResult:
    def __init__(self, tests, traits, covariates, candidates, h0, envs0, envs1):
        self._tests = tests
        self._traits = traits
        self._covariates = covariates
        self._candidates = candidates
        self._envs0 = envs0
        self._envs1 = envs1
        self._h0 = h0

    @property
    def stats(self):
        return self._dataframes["stats"].set_index("test")

    @property
    def effsizes(self):
        return self._dataframes["effsizes"]

    @property
    def _h0_dataframe(self):
        from pandas import DataFrame

        covariates = list(self._covariates)

        h0 = []
        for i, trait in enumerate(self._traits):
            for j, c in enumerate(covariates):
                h0.append([trait, "covariate", c, self._h0["effsizes"][i, j]])

        columns = ["trait", "effect_type", "effect_name", "effsize"]
        return DataFrame(h0, columns=columns)

    @property
    def _h1_dataframe(self):
        from pandas import DataFrame
        from itertools import product

        covariates = list(self._covariates)
        envs0 = list(self._envs0)

        h1 = []
        for i, test in enumerate(self._tests):
            candidates = list(self._candidates[test["idx"]])

            effsizes = test["h1"]["covariate_effsizes"]
            for j, trait in enumerate(self._traits):
                for l, c in enumerate(covariates):
                    v = [i, trait, "covariate", c, None, effsizes[j, l]]
                    h1.append(v)

            effsizes = test["h1"]["candidate_effsizes"]
            for j, trait in enumerate(self._traits):
                for l, ce in enumerate(product(candidates, envs0)):
                    env_name = "env0_" + str(ce[1])
                    v = [i, trait, "candidate", ce[0], env_name, effsizes[j, l]]
                    h1.append(v)

        columns = ["test", "trait", "effect_type", "effect_name", "env", "effsize"]
        return DataFrame(h1, columns=columns)

    @property
    def _h2_dataframe(self):
        from pandas import DataFrame
        from itertools import product

        envs0 = list(self._envs0)
        envs1 = list(self._envs1)
        covariates = list(self._covariates)

        h2 = []
        for i, test in enumerate(self._tests):
            candidates = list(self._candidates[test["idx"]])

            effsizes = test["h2"]["covariate_effsizes"]
            for j, trait in enumerate(self._traits):
                for l, c in enumerate(covariates):
                    v = [i, trait, "covariate", c, None, effsizes[j, l]]
                    h2.append(v)

            effsizes = test["h2"]["candidate_effsizes"]
            off = 0
            for j, trait in enumerate(self._traits):
                for l, ce in enumerate(product(candidates, envs0)):
                    env_name = "env0_" + str(ce[1])
                    v = [i, trait, "candidate", ce[0], env_name, effsizes[j, off + l]]
                    h2.append(v)

            off = len(candidates) * len(envs0)
            for j, trait in enumerate(self._traits):
                for l, ce in enumerate(product(candidates, envs1)):
                    env_name = "env1_" + str(ce[1])
                    v = [i, trait, "candidate", ce[0], env_name, effsizes[j, off + l]]
                    h2.append(v)

        columns = ["test", "trait", "effect_type", "effect_name", "env", "effsize"]
        return DataFrame(h2, columns=columns)

    @property
    def _stats_dataframe(self):
        from pandas import DataFrame

        stats = []
        for i, test in enumerate(self._tests):
            dof10 = test["h1"]["candidate_effsizes"].size
            dof20 = test["h2"]["candidate_effsizes"].size
            dof21 = dof20 - dof10
            stats.append(
                [
                    i,
                    self._h0["lml"],
                    test["h1"]["lml"],
                    test["h2"]["lml"],
                    dof10,
                    dof20,
                    dof21,
                ]
            )

        columns = ["test", "lml0", "lml1", "lml2", "dof10", "dof20", "dof21"]
        stats = DataFrame(stats, columns=columns)

        stats["pv10"] = lrt_pvalues(stats["lml0"], stats["lml1"], stats["dof10"])
        stats["pv20"] = lrt_pvalues(stats["lml0"], stats["lml2"], stats["dof20"])
        stats["pv21"] = lrt_pvalues(stats["lml1"], stats["lml2"], stats["dof21"])

        return stats

    @property
    @cache
    def _dataframes(self):
        h0 = self._h0_dataframe
        h1 = self._h1_dataframe
        h2 = self._h2_dataframe
        stats = self._stats_dataframe

        return {"stats": stats, "effsizes": {"h0": h0, "h1": h1, "h2": h2}}

    def __repr__(self):
        return ""
        # from textwrap import TextWrapper
        # from numpy import percentile, asarray

        # msg = "Null model\n"
        # msg += "----------\n\n"
        # v0 = self.background_variance
        # v1 = self.foreground_variance
        # msg += f"  y ~ 𝓝(M𝜶 + (E₀⊙Gᵢ)𝞫, {v0:.2f}*K + {v1:.2f}*I)\n"

        # prefix = "  M = "
        # s = " " * len(prefix)
        # wrapper = TextWrapper(initial_indent=prefix, width=88, subsequent_indent=s)
        # msg += wrapper.fill(str(self._covariates)) + "\n"

        # prefix = "  𝜶 = "
        # s = " " * len(prefix)
        # wrapper = TextWrapper(initial_indent=prefix, width=88, subsequent_indent=s)

        # effsizes = asarray(self.covariate_effsizes, float).ravel()
        # stats = self.stats

        # msg += wrapper.fill(str(effsizes)) + "\n"
        # msg += "  Log marg. lik.: {}\n".format(self._null_lml)
        # msg += "  Number of models: {}\n\n".format(stats.shape[0])
        # msg += "Alt model\n"
        # msg += "---------\n\n"
        # msg += f"  y ~ 𝓝(M𝜶 + (E₀⊙Gᵢ)𝞫₀ + (E₁⊙Gᵢ)𝞫₁, {v0:.2f}*K + {v1:.2f}*I)\n"

        # msg += "  Min. p-value: {}\n".format(min(stats["pvalue"]))
        # msg += "  First perc. p-value: {}\n".format(percentile(stats["pvalue"], 1))
        # msg += "  Max. log marg. lik.: {}\n".format(max(stats["alt lml"]))
        # msg += "  99th perc. log marg. lik.: {}\n".format(
        #     percentile(stats["alt lml"], 99)
        # )
        # msg += "  Number of models: {}".format(stats.shape[0])

        # return msg


# def _lik_formulae(lik):
#     msg = ""

#     if lik == "bernoulli":
#         msg += f"  yᵢ ~ Bern(μᵢ=g(zᵢ)), where g(x)=1/(1+e⁻ˣ)\n"
#     elif lik == "probit":
#         msg += f"  yᵢ ~ Bern(μᵢ=g(zᵢ)), where g(x)=Φ(x)\n"
#     elif lik == "binomial":
#         msg += f"  yᵢ ~ Binom(μᵢ=g(zᵢ), nᵢ), where g(x)=1/(1+e⁻ˣ)\n"
#     elif lik == "poisson":
#         msg += f"  yᵢ ~ Poisson(λᵢ=g(zᵢ)), where g(x)=eˣ\n"

#     return msg
