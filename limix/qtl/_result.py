from limix.stats import effsizes_se, lrt_pvalues
from functools import lru_cache


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
        msg += f"  y ~ 𝓝(M𝜶, {v0:.2f}*K + {v1:.2f}*I)\n"

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
        msg += f"  y ~ 𝓝(M𝜶 + Gᵢ, {v0:.2f}*K + {v1:.2f}*I)\n"

        msg += "  Min. p-value: {}\n".format(min(stats["pvalue"]))
        msg += "  First perc. p-value: {}\n".format(percentile(stats["pvalue"], 1))
        msg += "  Max. log marg. lik.: {}\n".format(max(stats["alt lml"]))
        msg += "  99th perc. log marg. lik.: {}\n".format(
            percentile(stats["alt lml"], 99)
        )
        msg += "  Number of models: {}\n".format(stats.shape[0])

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

                effsizes = effsizes[j + 1 :]
                for j, inter in enumerate(inters1):
                    alt.append([i, c, "inter1_" + str(inter), effsizes[j]])

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
        msg += "  Number of models: {}\n".format(stats.shape[0])

        return msg


# class QTLModel(object):
#     r"""Result of a QTL analysis.

#     An instance of this class is returned by :func:`limix.qtl.st_scan`.
#     """

#     def __init__(self, null_lml, alt_lmls, effsizes, null_covariate_effsizes):
#         self._null_lml = null_lml
#         self._alt_lmls = alt_lmls
#         self._effsizes = effsizes
#         self._null_covariate_effsizes = null_covariate_effsizes
#         alt_lmls.name = "alt lmls"
#         effsizes.name = "effsizes"

#     def _get_null_series(self):
#         from pandas import concat, Series

#         a = self._null_covariate_effsizes
#         b = Series(data=[self._null_lml], index=["null_lml"])
#         return concat([a, b])

#     def _get_alt_dataframe(self):
#         from pandas import DataFrame

#         df = DataFrame({"alt_lmls": self._alt_lmls, "effsizes": self._effsizes})
#         return df

#     @property
#     def null_lml(self):
#         r"""Log of the marginal likelihood under the null hypothesis.

#         Returns
#         -------
#         float
#             Log of marginal likelihood.
#         """
#         return self._null_lml

#     @property
#     def alt_lmls(self):
#         r"""Log of the marginal likelihoods across tested variants.

#         Returns
#         -------
#         array_like
#             Log of marginal likelihoods.
#         """
#         return self._alt_lmls

#     @property
#     def variant_effsizes(self):
#         r"""Variant effect-sizes.

#         Returns
#         -------
#         array_like
#             Estimated variant effect sizes.
#         """
#         return self._effsizes

#     @property
#     def variant_pvalues(self):
#         r"""Variant p-values.

#         Returns
#         -------
#         array_like
#             Association significance between variant and phenotype.
#         """
#         from xarray import zeros_like

#         pv = zeros_like(self._alt_lmls)
#         pv[:] = lrt_pvalues(self.null_lml, self._alt_lmls.values)
#         pv.name = "pv"

#         return pv

#     @property
#     def variant_effsizes_se(self):
#         r"""Standard errors of the variant effect sizes.

#         Returns
#         -------
#         array_like
#             Estimated standard errors of the variant effect sizes.
#         """
#         from xarray import zeros_like

#         ese = zeros_like(self._alt_lmls)
#         ese[:] = effsizes_se(
#             self.variant_effsizes.values.ravel(), self.variant_pvalues.values.ravel()
#         )
#         ese.name = "effsizes std"
#         return ese

#     @property
#     def null_covariate_effsizes(self):
#         r"""Covariate effect-sizes under the null hypothesis.

#         Returns
#         -------
#         array_like
#             Estimated covariant effect sizes under the null hypothesis.
#         """
#         return self._null_covariate_effsizes

#     def to_csv(self, path_or_buf_null, path_or_buf_alt):

#         null = self._get_null_series()
#         alt = self._get_alt_dataframe()

#         null.to_csv(path_or_buf_null, header=False)
#         alt.to_csv(path_or_buf_alt, header=False)

#     def __repr__(self):
#         import re
#         from pandas import DataFrame

#         data = dict(
#             effsizes=self.variant_effsizes.values.ravel(),
#             effsizes_se=self.variant_effsizes_se.values.ravel(),
#             pvalues=self.variant_pvalues.values.ravel(),
#         )

#         variant_msg = str(DataFrame(data=data).describe())

#         variant_lines = variant_msg.split("\n")

#         pline = variant_lines[1]
#         count_line = re.sub(r"(\d+)\.0", r" \1.", pline)
#         while pline != count_line:
#             pline = count_line
#             count_line = re.sub(r"(\d+)\.0", r" \1.", pline)

#         variant_lines[1] = re.sub(r"(\d+)\.", r" \1", count_line)
#         variant_msg = "\n".join(variant_lines)

#         data = self.null_covariate_effsizes
#         k = data.index.values
#         v = [[vi] for vi in data.values]

#         df = DataFrame(data=dict(zip(list(k), list(v)))).sort_index(axis=1)
#         covariate_msg = str(df)
#         covariate_msg = "\n".join([x[2:] for x in covariate_msg.split("\n")])

#         msg = "Variants\n--------\n" + variant_msg
#         msg += "\n\nCovariate effect sizes for H0\n"
#         msg += "-----------------------------\n"
#         msg += covariate_msg

#         return msg

# def __str__(self):
#     if PY2:
#         return self.__unicode__().encode("utf-8")
#     return self.__repr__()

# def __unicode__(self):
#     return self.__repr__()


# class STIScanResult(object):
#     def __init__(
#         self, null_lml, alt0_lmls, alt1_lmls, effsizes0, effsizes1, covariate_effsizes
#     ):
#         self._null_lml = null_lml
#         self._alt0_lmls = alt0_lmls
#         self._alt1_lmls = alt1_lmls
#         self._effsizes0 = effsizes0
#         self._effsizes1 = effsizes1
#         self._covariate_effsizes = covariate_effsizes
#         alt0_lmls.name = "alt0 lmls"
#         alt1_lmls.name = "alt1 lmls"
#         effsizes0.name = "effsizes0"
#         effsizes1.name = "effsizes1"

#     def _get_null_series(self):
#         from pandas import concat, Series

#         a = self._null_covariate_effsizes
#         b = Series(data=[self._null_lml], index=["null_lml"])
#         return concat([a, b])

#     def _get_alt_dataframe(self):
#         from pandas import DataFrame

#         df = DataFrame({"alt_lmls": self._alt_lmls, "effsizes": self._effsizes})
#         return df

#     @property
#     def null_lml(self):
#         r"""Log of the marginal likelihood under the null hypothesis.

#         Returns
#         -------
#         float
#             Log of marginal likelihood.
#         """
#         return self._null_lml

#     @property
#     def alt_lmls(self):
#         r"""Log of the marginal likelihoods across tested variants.

#         Returns
#         -------
#         array_like
#             Log of marginal likelihoods.
#         """
#         return self._alt_lmls

#     @property
#     def variant_effsizes(self):
#         r"""Variant effect-sizes.

#         Returns
#         -------
#         array_like
#             Estimated variant effect sizes.
#         """
#         return self._effsizes

#     @property
#     def variant_pvalues(self):
#         r"""Variant p-values.

#         Returns
#         -------
#         array_like
#             Association significance between variant and phenotype.
#         """
#         from xarray import zeros_like

#         pv = zeros_like(self._alt_lmls)
#         pv[:] = lrt_pvalues(self.null_lml, self._alt_lmls.values)
#         pv.name = "pv"

#         return pv

#     @property
#     def variant_effsizes_se(self):
#         r"""Standard errors of the variant effect sizes.

#         Returns
#         -------
#         array_like
#             Estimated standard errors of the variant effect sizes.
#         """
#         from xarray import zeros_like

#         ese = zeros_like(self._alt_lmls)
#         ese[:] = effsizes_se(
#             self.variant_effsizes.values.ravel(), self.variant_pvalues.values.ravel()
#         )
#         ese.name = "effsizes std"
#         return ese

#     @property
#     def null_covariate_effsizes(self):
#         r"""Covariate effect-sizes under the null hypothesis.

#         Returns
#         -------
#         array_like
#             Estimated covariant effect sizes under the null hypothesis.
#         """
#         return self._null_covariate_effsizes

#     def to_csv(self, path_or_buf_null, path_or_buf_alt):

#         null = self._get_null_series()
#         alt = self._get_alt_dataframe()

#         null.to_csv(path_or_buf_null, header=False)
#         alt.to_csv(path_or_buf_alt, header=False)

#     def __repr__(self):
#         import re
#         from pandas import DataFrame

#         data = dict(
#             effsizes=self.variant_effsizes.values.ravel(),
#             effsizes_se=self.variant_effsizes_se.values.ravel(),
#             pvalues=self.variant_pvalues.values.ravel(),
#         )

#         variant_msg = str(DataFrame(data=data).describe())

#         variant_lines = variant_msg.split("\n")

#         pline = variant_lines[1]
#         count_line = re.sub(r"(\d+)\.0", r" \1.", pline)
#         while pline != count_line:
#             pline = count_line
#             count_line = re.sub(r"(\d+)\.0", r" \1.", pline)

#         variant_lines[1] = re.sub(r"(\d+)\.", r" \1", count_line)
#         variant_msg = "\n".join(variant_lines)

#         data = self.null_covariate_effsizes
#         k = data.index.values
#         v = [[vi] for vi in data.values]

#         df = DataFrame(data=dict(zip(list(k), list(v)))).sort_index(axis=1)
#         covariate_msg = str(df)
#         covariate_msg = "\n".join([x[2:] for x in covariate_msg.split("\n")])

#         msg = "Variants\n--------\n" + variant_msg
#         msg += "\n\nCovariate effect sizes for H0\n"
#         msg += "-----------------------------\n"
#         msg += covariate_msg

#         return msg

#     def __str__(self):
#         if py2:
#             return self.__unicode__().encode("utf-8")
#         return self.__repr__()

#     def __unicode__(self):
#         return self.__repr__()


def _make_sure_iterabble(x):
    from collections.abc import Iterable

    if not isinstance(x, Iterable):
        return [x]

    return x


if __name__ == "__main__":
    import numpy as np
    from limix._data import conform_dataset
    from limix._data import asarray as _asarray

    random = np.random.RandomState(0)
    G = random.randn(5, 5)
    E0 = random.randn(5, 1)
    E1 = random.randn(5, 2)
    y = random.randn(5)
    M = random.randn(5, 2)
    K = random.randn(5, 5)
    K = K @ K.T

    data = conform_dataset(y, M, G, K)
    y = data["y"]
    M = data["M"]
    M = M.assign_coords(covariate=["offset", 0])
    G = data["G"]
    K = data["K"]
    E0 = _asarray(E0, "inter0", ["sample", "inter"])
    E1 = _asarray(E1, "inter1", ["sample", "inter"])

    # r = ST_ScanResultFactory("normal", M.covariate, G.candidate)

    # r.set_null(-10.203, random.randn(2), 0.2, 0.6)
    # r.add_test(0, random.randn(1), -10.0)
    # r.add_test([1, 2], random.randn(2), -9.0)

    # r = r.create()

    # print(r.foreground_variance)
    # print(r.background_variance)
    # print(r.covariate_effsizes)
    # print()
    # print(r.alt_effsizes)
    # print()
    # print(r.stats)
    # print()
    # print(r)

    r = ST_IScanResultFactory(M.covariate, G.candidate, E1.inter, E0.inter)
    r.set_null(-10.203, random.randn(2), 0.2, 0.6)
    r.add_test(slice(0, 3), random.randn(1 * 3), -10.2, random.randn(3 * 3), -9.0)
    r.add_test([3, 4], random.randn(1 * 2), -10.1, random.randn(3 * 2), -9.2)

    r = r.create()

    print(r.foreground_variance)
    print(r.background_variance)
    print()
    print(r.covariate_effsizes)
    print()
    print(r.null_effsizes)
    print()
    print(r.alt_effsizes)
    print()
    print(r.stats)
    print()
    print(r)
