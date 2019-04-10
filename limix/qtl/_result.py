from limix._bits import unvec
from limix._cache import cache
from limix.stats import lrt_pvalues

from ._aligned import Aligned
from ._table import Table


class SModelResult:
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


class ScanResultFactory:
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
                "covariate_effsizes": _2d_shape(h["covariate_effsizes"]),
                "candidate_effsizes": _2d_shape(h["candidate_effsizes"]),
                "covariate_effsizes_se": _2d_shape(h["covariate_effsizes_se"]),
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

        h1 = []
        for i, test in enumerate(self._tests):
            candidates = list(self._candidates[test["idx"]])

            effsizes = test["h1"]["covariate_effsizes"]
            effsizes_se = test["h1"]["covariate_effsizes_se"]
            for j, trait in enumerate(self._traits):
                for l, c in enumerate(covariates):
                    eff = effsizes[l, j]
                    eff_se = effsizes_se[l, j]
                    v = [i, str(trait), "covariate", str(c), None, eff, eff_se]
                    h1.append(v)

            effsizes = test["h1"]["candidate_effsizes"]
            effsizes_se = test["h1"]["candidate_effsizes_se"]
            for j, e in enumerate(envs0):
                for l, c in enumerate(candidates):
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

        h2 = []
        for i, test in enumerate(self._tests):
            candidates = list(self._candidates[test["idx"]])

            effsizes = test["h2"]["covariate_effsizes"]
            effsizes_se = test["h2"]["covariate_effsizes_se"]
            for j, trait in enumerate(self._traits):
                for l, c in enumerate(covariates):
                    eff = effsizes[l, j]
                    eff_se = effsizes_se[l, j]
                    v = [i, str(trait), "covariate", str(c), None, eff, eff_se]
                    h2.append(v)

            effsizes = test["h2"]["candidate_effsizes"]
            effsizes_se = test["h2"]["candidate_effsizes_se"]
            off = 0
            for j, e in enumerate(envs0):
                for l, c in enumerate(candidates):
                    env_name = "env0_" + str(e)
                    eff = effsizes[l, off + j]
                    eff_se = effsizes_se[l, off + j]
                    v = [i, str(trait), "candidate", str(c), env_name, eff, eff_se]
                    h2.append(v)

            off = len(envs0)
            for j, e in enumerate(envs1):
                for l, c in enumerate(candidates):
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

    def _repr_three_hypothesis(self):
        from numpy import asarray

        lik = self._h0.likelihood
        covariates = self._covariates
        lml = self._h0.lml
        effsizes = asarray(self.h0.effsizes["effsize"], float).ravel()
        effsizes_se = asarray(self.h0.effsizes["effsize_se"], float).ravel()
        stats = self.stats
        # v0 = self.h0.variances["fore_covariance"]
        # v1 = self.h0.variances["back_covariance"]

        msg = draw_title("Hypothesis 0")
        msg += draw_model(lik, "(A⊗𝙼)𝜶", "C₀⊗𝙺 + C₁⊗𝙸")
        msg += _draw_hyp0_summary(covariates, effsizes, effsizes_se, lml)

        msg += draw_title("Hypothesis 1")
        msg += draw_model(lik, "(A⊗𝙼)𝜶 + (A₀⊗G)𝜶₀", "s(C₀⊗𝙺 + C₁⊗𝙸)")
        msg += draw_alt_hyp_table(1, self.stats, self.effsizes)

        msg += draw_title("Hypothesis 2")
        msg += draw_model(lik, "(A⊗𝙼)𝜶 + (A₀⊗G)𝜶₀ + (A₁⊗G)𝜶₁", "s(C₀⊗𝙺 + C₁⊗𝙸)")
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
        # v0 = self.h0.variances["fore_covariance"].item()
        # v1 = self.h0.variances["back_covariance"].item()

        covariance = self._covariance_expr()

        msg = draw_title("Hypothesis 0")
        msg += draw_model(lik, "(A⊗𝙼)𝜶", "C₀⊗𝙺 + C₁⊗𝙸")
        msg += _draw_hyp0_summary(covariates, effsizes, effsizes_se, lml)

        if alt_hyp == 1:
            mean = "𝙼𝜶 + (A₀⊗G)𝜶₀"
            col = "𝓗₀ vs 𝓗₁"
        else:
            mean = "𝙼𝜶 + (A₁⊗G)𝜶₁"
            col = "𝓗₀ vs 𝓗₂"

        msg += draw_title(f"Hypothesis {alt_hyp}")
        msg += draw_model(lik, mean, f"s({covariance})")
        msg += draw_alt_hyp_table(alt_hyp, self.stats, self.effsizes)

        msg += draw_title("Likelihood-ratio test p-values")
        msg += draw_lrt_table([col], [f"pv{alt_hyp}0"], stats)
        return msg

    def __repr__(self):
        if len(self._envs0) > 0 and len(self._envs0) == len(self._envs1):
            return self._repr_two_hypothesis(1)
        elif len(self._envs0) > 0 and len(self._envs1) > 0:
            return self._repr_three_hypothesis()
        elif len(self._envs0) == 0 and len(self._envs1) > len(self._envs0):
            return self._repr_two_hypothesis(2)
        raise ValueError("There is no alternative hypothesis.")


def draw_title(title):
    msg = f"{title}\n"
    msg += "=" * len(title) + "\n\n"
    return msg


def draw_model(lik, mean, covariance):
    if lik == "normal":
        var = "𝐲"
    else:
        var = "𝐳"

    msg = f"{var} ~ 𝓝({mean}, {covariance})"
    msg += _lik_formulae(lik)
    return msg


def draw_alt_hyp_table(hyp_num, stats, effsizes):
    cols = ["lml", "cov. effsizes", "cand. effsizes"]
    table = Table(cols, index=_describe_index())
    table.add_column(_describe(stats, f"lml{hyp_num}"))
    df = effsizes[f"h{hyp_num}"]
    table.add_column(_describe(df[df["effect_type"] == "covariate"], "effsize"))
    table.add_column(_describe(df[df["effect_type"] == "candidate"], "effsize"))
    return "\n" + table.draw() + "\n\n"


def draw_lrt_table(test_titles, pv_names, stats):
    table = Table(test_titles, index=_describe_index())

    for name in pv_names:
        pv = stats[name].describe().iloc[1:]
        table.add_column(pv)

    return table.draw()


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


def _describe(df, field):
    return df[field].describe().iloc[1:]


def _describe_index():
    return ["mean", "std", "min", "25%", "50%", "75%", "max"]


def _draw_hyp0_summary(covariates, effsizes, effsizes_se, lml):
    aligned = Aligned()
    aligned.add_item("M", covariates)
    aligned.add_item("𝜶", effsizes)
    aligned.add_item("se(𝜶)", effsizes_se)
    aligned.add_item("lml", lml)
    return aligned.draw() + "\n"
