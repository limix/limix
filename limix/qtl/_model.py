from __future__ import division

import sys

from pandas import DataFrame, Series

from limix.stats import effsizes_se, lrt_pvalues

if sys.version_info < (3, 0):
    PY2 = True
else:
    PY2 = False


class QTLModel(object):
    r"""Result of a QTL analysis.

    An instance of this class is returned by :func:`limix.qtl.scan`.
    """

    def __init__(self, null_lml, alt_lmls, effsizes, null_covariate_effsizes):
        self._null_lml = null_lml
        self._alt_lmls = alt_lmls
        self._effsizes = effsizes
        self._null_covariate_effsizes = null_covariate_effsizes

    @property
    def null_lml(self):
        r"""Log of the marginal likelihood under the null hypothesis.

        Returns
        -------
        float
            Log of marginal likelihood.
        """
        return self._null_lml

    @property
    def alt_lmls(self):
        r"""Log of the marginal likelihoods across tested variants.

        Returns
        -------
        array_like
            Log of marginal likelihoods.
        """
        return self._alt_lmls

    @property
    def variant_effsizes(self):
        r"""Variant effect-sizes.

        Returns
        -------
        array_like
            Estimated variant effect sizes.
        """
        return self._effsizes

    @property
    def variant_pvalues(self):
        r"""Variant p-values.

        Returns
        -------
        array_like
            Association significance between variant and phenotype.
        """
        pv = lrt_pvalues(self.null_lml, self.alt_lmls.values)
        return Series(pv, list(self.alt_lmls.keys()))

    @property
    def variant_effsizes_se(self):
        r"""Standard errors of the variant effect sizes.

        Returns
        -------
        array_like
            Estimated standard errors of the variant effect sizes.
        """
        ese = effsizes_se(self.variant_effsizes.values, self.variant_pvalues.values)
        return Series(ese, list(self.alt_lmls.keys()))

    @property
    def null_covariate_effsizes(self):
        r"""Covariate effect-sizes under the null hypothesis.

        Returns
        -------
        array_like
            Estimated covariant effect sizes under the null hypothesis.
        """
        return self._null_covariate_effsizes

    def __repr__(self):
        data = dict(
            effsizes=self.variant_effsizes,
            effsizes_se=self.variant_effsizes_se,
            pvalues=self.variant_pvalues,
        )

        variant_msg = str(DataFrame(data=data).describe())

        data = self.null_covariate_effsizes
        k = data.index.values
        v = [[vi] for vi in data.values]

        df = DataFrame(data=dict(zip(list(k), list(v)))).sort_index(axis=1)
        covariate_msg = str(df)
        covariate_msg = "\n".join([x[2:] for x in covariate_msg.split("\n")])

        msg = "Variants\n--------\n" + variant_msg
        msg += "\n\nCovariate effect sizes for H0\n"
        msg += "-----------------------------\n"
        msg += covariate_msg

        return msg

    def __str__(self):
        if PY2:
            return self.__unicode__().encode("utf-8")
        return self.__repr__()

    def __unicode__(self):
        return self.__repr__()
