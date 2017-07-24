from __future__ import division

from numpy import abs as npy_abs
from numpy import sqrt
from scipy.stats import chi2
from tabulate import tabulate
from pandas import Series

from limix.stats import effsizes_se, lrt_pvalues


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
        ese = effsizes_se(self.variant_effsizes.values,
                          self.variant_pvalues.values)
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

    def __str__(self):
        from pandas import DataFrame

        data = dict(
            effsizes=self.variant_effsizes,
            effsizes_se=self.variant_effsizes_se,
            pvalues=self.variant_pvalues)

        variant_msg = str(DataFrame(data=data).describe())

        data = self.null_covariate_effsizes

        covariate_msg = tabulate(
            [data.values], headers=list(data.keys()), tablefmt="plain")

        msg = 'Variants\n' + variant_msg
        msg += '\n\nCovariate effect sizes for the'
        msg += ' null model\n' + covariate_msg

        return msg
