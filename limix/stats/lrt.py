from __future__ import division

from numpy import abs as npy_abs
from numpy import asarray, sqrt
from scipy.stats import chi2


def lrt_pvalues(null_lml, alt_lmls, dof=1):
    r"""Compute p-values from likelihood ratios.

    These are likelihood ratio test p-values.

    Parameters
    ----------
    null_lml : float
        Log of the marginal likelihood under the null hypothesis.
    alt_lmls : array_like
        Log of the marginal likelihoods under the alternative hypotheses.
    dof : int
        Degrees of freedom.

    Returns
    -------
    array_like
        P-values.
    """
    from scipy.stats import chi2
    lrs = -2 * null_lml + 2 * asarray(alt_lmls)
    return chi2(df=dof).sf(lrs)


def effsizes_se(effsizes, pvalues):
    r"""Standard errors of the effect sizes.

    Parameters
    ----------
    effsizes : array_like
        Effect sizes.
    pvalues : array_like
        Association significance corresponding to those effect sizes.

    Returns
    -------
    array_like
        Standard errors of the effect sizes.
    """
    return npy_abs(effsizes) / sqrt(chi2(1).isf(pvalues))
