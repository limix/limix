from __future__ import division


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
    from numpy_sugar import epsilon
    from numpy import asarray, clip, inf

    lrs = clip(-2 * null_lml + 2 * asarray(alt_lmls), epsilon.super_tiny, inf)
    pv = chi2(df=dof).sf(lrs)
    return clip(pv, epsilon.super_tiny, 1 - epsilon.tiny)


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
    from scipy.stats import chi2
    from numpy import abs as npy_abs
    from numpy import sqrt

    return npy_abs(effsizes) / sqrt(chi2(1).isf(pvalues))
