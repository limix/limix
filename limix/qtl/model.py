from __future__ import division

from numpy import abs as npy_abs
from numpy import sqrt
from scipy.stats import chi2
from tabulate import tabulate

from limix.stats import effsizes_se, lrt_pvalues


class QTLModel(object):
    def __init__(self):
        pass

    @property
    def null_lml(self):
        raise NotImplementedError

    @property
    def alt_lmls(self):
        raise NotImplementedError

    @property
    def variant_effsizes(self):
        r"""
        Returns
        -------
        array_like
            Estimated variant effect sizes.
        """
        raise NotImplementedError

    @property
    def variant_pvalues(self):
        r"""
        Returns
        -------
            Association significance between variant and phenotype.
        """
        raise NotImplementedError

    @property
    def variant_effsizes_se(self):
        r"""
        Returns
        -------
        array_like
            Estimated standard errors of the variant effect sizes.
        """
        effsizes = self.variant_effsizes
        pv = self.variant_pvalues
        return npy_abs(effsizes) / sqrt(chi2(1).isf(pv))

    @property
    def null_covariant_effsizes(self):
        r"""
        Returns
        -------
        array_like
            Estimated covariant effect sizes under the null hypothesis.
        """
        raise NotImplementedError


class QTLModel_LMM(QTLModel):
    def __init__(self, null_lml, alt_lmls, effsizes, null_covariate_effsizes):
        self._null_lml = null_lml
        self._alt_lmls = alt_lmls
        self._effsizes = effsizes
        self._null_covariate_effsizes = null_covariate_effsizes

    @property
    def null_lml(self):
        return self._null_lml

    @property
    def alt_lmls(self):
        return self._alt_lmls

    @property
    def variant_effsizes(self):
        return self._effsizes

    @property
    def variant_effsizes_se(self):
        return effsizes_se(self.variant_effsizes, self.variant_pvalues)

    @property
    def variant_pvalues(self):
        return lrt_pvalues(self.null_lml, self.alt_lmls)

    @property
    def null_covariate_effsizes(self):
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
            [list(data.values())], headers=list(data.keys()), tablefmt="plain")

        msg = 'Variants\n' + variant_msg
        msg += '\n\nCovariate effect sizes for the'
        msg += ' null model\n' + covariate_msg

        return msg


class QTLModel_GLMM(QTLModel):
    def __init__(self, null_lml, alt_lmls, effsizes, null_covariate_effsizes):
        self._null_lml = null_lml
        self._alt_lmls = alt_lmls
        self._effsizes = effsizes
        self._null_covariate_effsizes = null_covariate_effsizes

    @property
    def null_lml(self):
        return self._null_lml

    @property
    def alt_lmls(self):
        return self._alt_lmls

    @property
    def variant_effsizes(self):
        return self._effsizes

    @property
    def variant_effsizes_se(self):
        return effsizes_se(self.variant_effsizes, self.variant_pvalues)

    @property
    def variant_pvalues(self):
        return lrt_pvalues(self.null_lml, self.alt_lmls)

    @property
    def null_covariate_effsizes(self):
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
            [list(data.values())], headers=list(data.keys()), tablefmt="plain")

        msg = 'Variants\n' + variant_msg
        msg += '\n\nCovariate effect sizes for the'
        msg += ' null model\n' + covariate_msg

        return msg
