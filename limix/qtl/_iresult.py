from limix._cache import cache
from limix.stats import lrt_pvalues

from ._result import SModelResult
from ._table import Table


class IScanResultFactory:
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

        self._h0 = SModelResult(lik, traits, covariates, lml, beta, beta_se, C0, C1)
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
            return {
                "lml": float(h["lml"]),
                "covariate_effsizes": asarray(h["covariate_effsizes"], float),
                "candidate_effsizes": _2d_shape(h["candidate_effsizes"]),
                "covariate_effsizes_se": asarray(h["covariate_effsizes_se"], float),
                "candidate_effsizes_se": _2d_shape(h["candidate_effsizes_se"]),
                "scale": float(h["scale"]),
            }

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
        """
        TODO
        """
        return self._dataframes["stats"].set_index("test")

    @property
    def effsizes(self):
        """
        TODO
        """
        return self._dataframes["effsizes"]

    @property
    def h0(self):
        """
        TODO
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

    def _repr_two_hypothesis(self):
        from numpy import asarray, isnan

        lik = self._h0.likelihood
        effsizes = asarray(self.h0.effsizes["effsize"], float).ravel()
        effsizes_se = asarray(self.h0.effsizes["effsize_se"], float).ravel()
        stats = self.stats

        msg = "Null model\n"
        msg += "----------\n\n"
        v0 = self.h0.variances["fore_covariance"].item()
        v1 = self.h0.variances["back_covariance"].item()
        if lik == "normal":
            var = "𝐲"
        else:
            var = "𝐳"
        if isnan(v0):
            msg += f"{var} ~ 𝓝(𝙼𝜶, {v1:.4f}⋅𝙸)"
        else:
            msg += f"{var} ~ 𝓝(𝙼𝜶, {v0:.4f}⋅𝙺 + {v1:.4f}⋅𝙸)"
        msg += _lik_formulae(self._h0.likelihood)

        msg += _item_repr("𝙼     = ", self._covariates)
        msg += _item_repr("𝜶     = ", effsizes)
        msg += _item_repr("se(𝜶) = ", effsizes_se)

        msg += "lml   = {}\n\n".format(self.h0.lml)

        msg += "Alt model\n"
        msg += "---------\n\n"
        if isnan(v0):
            msg += f"{var} ~ 𝓝(𝙼𝜶 + 𝙶𝞫, {v1:.4f}⋅𝙸)"
        else:
            msg += f"{var} ~ 𝓝(𝙼𝜶 + 𝙶𝞫, {v0:.4f}⋅𝙺 + {v1:.4f}⋅𝙸)"
        msg += _lik_formulae(self._h0.likelihood)

        msg += "min(pv)  = {}\n".format(min(stats["pv20"]))
        msg += "max(lml) = {}\n".format(max(stats["lml2"]))

        return msg

    def _repr_three_hypothesis(self):
        from numpy import asarray

        lik = self._h0.likelihood
        covariates = self._covariates
        lml = self._h0.lml
        effsizes = asarray(self.h0.effsizes["effsize"], float).ravel()
        effsizes_se = asarray(self.h0.effsizes["effsize_se"], float).ravel()
        stats = self.stats
        v0 = self.h0.variances["fore_covariance"].item()
        v1 = self.h0.variances["back_covariance"].item()

        msg = _draw_hypothesis_zero(lik, v0, v1, covariates, effsizes, effsizes_se, lml)

        msg += _section("Hypothesis 1", " + (𝙶⊙𝙴₀)𝛃₀", lik, v0, v1, True)
        msg += _draw_alt_hypothesis_table(1, self.stats, self.effsizes)

        msg += _section("Hypothesis 2", " + (𝙶⊙𝙴₀)𝛃₀ + (𝙶⊙𝙴₁)𝛃₁", lik, v0, v1, True)
        msg += _draw_alt_hypothesis_table(2, self.stats, self.effsizes)

        msg += _draw_lrt_section(
            ["𝓗₀ vs 𝓗₁", "𝓗₀ vs 𝓗₂", "𝓗₁ vs 𝓗₂"], ["pv10", "pv20", "pv21"], stats
        )
        return msg

    def __repr__(self):
        if len(self._envs0) > 0 and len(self._envs1) == 0:
            return self._repr_two_hypothesis()
        elif len(self._envs0) > 0 and len(self._envs1) > 0:
            return self._repr_three_hypothesis()
        elif len(self._envs0) == 0 and len(self._envs1) > 0:
            return self._repr_two_hypothesis()
        raise ValueError("There is no environment to interact with.")


def _item_repr(prefix, item):
    from textwrap import TextWrapper

    s = " " * len(prefix)
    wrapper = TextWrapper(initial_indent=prefix, width=88, subsequent_indent=s)
    return wrapper.fill(str(item)) + "\n"


def _lik_formulae(lik):
    msg = ""

    if lik == "bernoulli":
        msg += f" for yᵢ ~ Bern(μᵢ=g(zᵢ)) and g(x)=1/(1+e⁻ˣ)\n"
    elif lik == "probit":
        msg += f" for yᵢ ~ Bern(μᵢ=g(zᵢ)) and g(x)=Φ(x)\n"
    elif lik == "binomial":
        msg += f" for yᵢ ~ Binom(μᵢ=g(zᵢ), nᵢ) and g(x)=1/(1+e⁻ˣ)\n"
    elif lik == "poisson":
        msg += f" for yᵢ ~ Poisson(λᵢ=g(zᵢ)) and g(x)=eˣ\n"
    else:
        msg += "\n"
    return msg


def _title(title):
    msg = f"{title}\n"
    msg += "=" * len(title) + "\n\n"
    return msg


def _section(model_name, cand_term, lik, v0, v1, scale):
    from numpy import isnan

    msg = _title(model_name)
    if lik == "normal":
        var = "𝐲"
    else:
        var = "𝐳"

    if cand_term is None:
        cand_term = ""

    if scale:
        left = "s("
        right = ")"
    else:
        left = ""
        right = ""

    if isnan(v0):
        msg += f"{var} ~ 𝓝(𝙼𝜶{cand_term}, {left}{v1:.4f}⋅𝙸{right})"
    else:
        msg += f"{var} ~ 𝓝(𝙼𝜶{cand_term}, {left}{v0:.4f}⋅𝙺 + {v1:.4f}⋅𝙸{right})"
    msg += _lik_formulae(lik)
    return msg


def _draw_hypothesis_zero(lik, v0, v1, covariates, effsizes, effsizes_se, lml):
    from ._aligned import Aligned

    msg = _section("Hypothesis 0", None, lik, v0, v1, False)
    aligned = Aligned()
    aligned.add_item("M", covariates)
    aligned.add_item("𝜶", effsizes)
    aligned.add_item("se(𝜶)", effsizes_se)
    aligned.add_item("lml", lml)
    msg += aligned.draw() + "\n"
    return msg


def _describe(df, field):
    return df[field].describe().iloc[1:]


def _describe_index():
    return ["mean", "std", "min", "25%", "50%", "75%", "max"]


def _draw_alt_hypothesis_table(hyp_num, stats, effsizes):
    cols = ["lml", "cov. effsizes", "cand. effsizes"]
    table = Table(cols, index=_describe_index())
    table.add_column(_describe(stats, f"lml{hyp_num}"))
    df = effsizes[f"h{hyp_num}"]
    table.add_column(_describe(df[df["effect_type"] == "covariate"], "effsize"))
    table.add_column(_describe(df[df["effect_type"] == "candidate"], "effsize"))
    return "\n" + table.draw() + "\n\n"


def _draw_lrt_section(test_titles, pv_names, stats):
    msg = _title("Likelihood-ratio test p-values")

    table = Table(test_titles, index=_describe_index())

    for name in pv_names:
        pv = stats[name].describe().iloc[1:]
        table.add_column(pv)

    msg += table.draw()
    return msg
