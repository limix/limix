from limix._cache import cache
from limix.stats import lrt_pvalues

from ._aligned import Aligned
from ._draw import draw_alt_hyp_table, draw_lrt_table, draw_model, draw_title


class IScanResult:
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
        """
        Statistics.
        """
        return self._dataframes["stats"].set_index("test")

    @property
    def effsizes(self):
        """
        Effect sizes.
        """
        return self._dataframes["effsizes"]

    @property
    def h0(self):
        """
        Hypothesis zero.
        """
        return self._h0

    @property
    def _h0_dataframe(self):
        from pandas import DataFrame

        covariates = list(self._covariates)

        h0 = []
        for i, trait in enumerate(self._traits):
            for j, c in enumerate(covariates):
                eff = self._h0["effsizes"][j, i]
                eff_se = self._h0["effsizes_se"][j, i]
                h0.append([trait, "covariate", c, eff, eff_se])

        columns = ["trait", "effect_type", "effect_name", "effsize", "effsize_se"]
        return DataFrame(h0, columns=columns)

    @property
    def _h1_dataframe(self):
        from pandas import DataFrame

        covariates = list(self._covariates)
        envs0 = list(self._envs0)

        trait = self._traits[0]
        h1 = []
        for i, test in enumerate(self._tests):
            candidates = list(self._candidates[test["idx"]])

            effsizes = test["h1"]["covariate_effsizes"]
            effsizes_se = test["h1"]["covariate_effsizes_se"]
            for j, c in enumerate(covariates):
                eff = effsizes[j]
                eff_se = effsizes_se[j]
                v = [i, str(trait), "covariate", str(c), None, eff, eff_se]
                h1.append(v)

            effsizes = test["h1"]["candidate_effsizes"]
            effsizes_se = test["h1"]["candidate_effsizes_se"]
            for l, c in enumerate(candidates):
                for j, e in enumerate(envs0):
                    env_name = "env0_" + str(e)
                    eff = effsizes[l, j]
                    eff_se = effsizes_se[l, j]
                    v = [i, str(trait), "candidate", str(c), env_name, eff, eff_se]
                    h1.append(v)

        columns = [
            "test",
            "trait",
            "effect_type",
            "effect_name",
            "env",
            "effsize",
            "effsize_se",
        ]
        return DataFrame(h1, columns=columns)

    @property
    def _h2_dataframe(self):
        from pandas import DataFrame

        envs0 = list(self._envs0)
        envs1 = list(self._envs1)
        covariates = list(self._covariates)

        trait = self._traits[0]
        h2 = []
        for i, test in enumerate(self._tests):
            candidates = list(self._candidates[test["idx"]])

            effsizes = test["h2"]["covariate_effsizes"]
            effsizes_se = test["h2"]["covariate_effsizes_se"]
            for j, c in enumerate(covariates):
                eff = effsizes[j]
                eff_se = effsizes_se[j]
                v = [i, str(trait), "covariate", str(c), None, eff, eff_se]
                h2.append(v)

            effsizes = test["h2"]["candidate_effsizes"]
            effsizes_se = test["h2"]["candidate_effsizes_se"]
            off = 0
            for l, c in enumerate(candidates):
                for j, e in enumerate(envs0):
                    env_name = "env0_" + str(e)
                    eff = effsizes[l, off + j]
                    eff_se = effsizes_se[l, off + j]
                    v = [i, str(trait), "candidate", str(c), env_name, eff, eff_se]
                    h2.append(v)

            off = len(envs0)
            for l, c in enumerate(candidates):
                for j, e in enumerate(envs1):
                    env_name = "env1_" + str(e)
                    eff = effsizes[l, off + j]
                    eff_se = effsizes_se[l, off + j]
                    v = [i, str(trait), "candidate", str(c), env_name, eff, eff_se]
                    h2.append(v)

        columns = [
            "test",
            "trait",
            "effect_type",
            "effect_name",
            "env",
            "effsize",
            "effsize_se",
        ]
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
                    self._h0.lml,
                    test["h1"]["lml"],
                    test["h2"]["lml"],
                    dof10,
                    dof20,
                    dof21,
                    test["h1"]["scale"],
                    test["h2"]["scale"],
                ]
            )

        columns = [
            "test",
            "lml0",
            "lml1",
            "lml2",
            "dof10",
            "dof20",
            "dof21",
            "scale1",
            "scale2",
        ]
        stats = DataFrame(stats, columns=columns)

        stats["pv10"] = lrt_pvalues(stats["lml0"], stats["lml1"], stats["dof10"])
        stats["pv20"] = lrt_pvalues(stats["lml0"], stats["lml2"], stats["dof20"])
        stats["pv21"] = lrt_pvalues(stats["lml1"], stats["lml2"], stats["dof21"])

        return stats

    @property
    @cache
    def _dataframes(self):
        h1 = self._h1_dataframe
        h2 = self._h2_dataframe
        stats = self._stats_dataframe

        return {"stats": stats, "effsizes": {"h1": h1, "h2": h2}}

    def _covariance_expr(self):
        from numpy import isnan

        v0 = self.h0.variances["fore_covariance"].item()
        v1 = self.h0.variances["back_covariance"].item()

        if isnan(v0):
            covariance = f"{v1:.3f}⋅𝙸"
        else:
            covariance = f"{v0:.3f}⋅𝙺 + {v1:.3f}⋅𝙸"

        return covariance

    def _repr_three_hypothesis(self):
        from numpy import asarray

        lik = self._h0.likelihood
        covariates = self._covariates
        lml = self._h0.lml
        effsizes = asarray(self.h0.effsizes["effsize"], float).ravel()
        effsizes_se = asarray(self.h0.effsizes["effsize_se"], float).ravel()
        stats = self.stats

        covariance = self._covariance_expr()

        msg = draw_title("Hypothesis 0")
        msg += draw_model(lik, "𝙼𝜶", covariance) + "\n"
        msg += _draw_hyp0_summary(covariates, effsizes, effsizes_se, lml)

        msg += draw_title("Hypothesis 1")
        msg += draw_model(lik, "𝙼𝜶 + (𝙶⊙𝙴₀)𝛃₀", f"s({covariance})")
        msg += draw_alt_hyp_table(1, self.stats, self.effsizes)

        msg += draw_title("Hypothesis 1")
        msg += draw_model(lik, "𝙼𝜶 + (𝙶⊙𝙴₀)𝛃₀ + (𝙶⊙𝙴₁)𝛃₁", f"s({covariance})")
        msg += draw_alt_hyp_table(2, self.stats, self.effsizes)

        msg += draw_title("Likelihood-ratio test p-values")
        cols = ["𝓗₀ vs 𝓗₁", "𝓗₀ vs 𝓗₂", "𝓗₁ vs 𝓗₂"]
        msg += draw_lrt_table(cols, ["pv10", "pv20", "pv21"], stats)
        return msg

    def _repr_two_hypothesis(self, alt_hyp):
        from numpy import asarray

        lik = self._h0.likelihood
        covariates = self._covariates
        lml = self._h0.lml
        effsizes = asarray(self.h0.effsizes["effsize"], float).ravel()
        effsizes_se = asarray(self.h0.effsizes["effsize_se"], float).ravel()
        stats = self.stats

        covariance = self._covariance_expr()

        msg = draw_title("Hypothesis 0")
        msg += draw_model(lik, "𝙼𝜶", covariance) + "\n"
        msg += _draw_hyp0_summary(covariates, effsizes, effsizes_se, lml)

        if alt_hyp == 1:
            mean = "𝙼𝜶 + (𝙶⊙𝙴₀)𝛃₀"
            col = "𝓗₀ vs 𝓗₁"
        else:
            mean = "𝙼𝜶 + (𝙶⊙𝙴₁)𝛃₁"
            col = "𝓗₀ vs 𝓗₂"

        msg += draw_title(f"Hypothesis {alt_hyp}")
        msg += draw_model(lik, mean, f"s({covariance})")
        msg += draw_alt_hyp_table(alt_hyp, self.stats, self.effsizes)

        msg += draw_title("Likelihood-ratio test p-values")
        msg += draw_lrt_table([col], [f"pv{alt_hyp}0"], stats)
        return msg

    def __repr__(self):
        if len(self._envs0) == 0:
            return self._repr_two_hypothesis(2)
        return self._repr_three_hypothesis()


def _draw_hyp0_summary(covariates, effsizes, effsizes_se, lml):
    aligned = Aligned()
    aligned.add_item("M", covariates)
    aligned.add_item("𝜶", effsizes)
    aligned.add_item("se(𝜶)", effsizes_se)
    aligned.add_item("lml", lml)
    return aligned.draw() + "\n"
