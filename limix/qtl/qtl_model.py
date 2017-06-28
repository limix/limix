from __future__ import division

from numpy import abs as npy_abs
from numpy import sqrt
from scipy.stats import chi2


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

    @property
    def alt_covariant_effsizes(self):
        r"""
        Returns
        -------
        array_like
            Estimated covariant effect sizes under the alternative hypothesis.
        """
        raise NotImplementedError
