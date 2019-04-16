from __future__ import division


def lrt_pvalues(null_lml, alt_lmls, dof=1):
    """
    Compute p-values from likelihood ratios.

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
    pvalues : ndarray
        P-values.
    """
    from scipy.stats import chi2
    from numpy_sugar import epsilon
    from numpy import asarray, clip, inf

    lrs = clip(-2 * null_lml + 2 * asarray(alt_lmls, float), epsilon.super_tiny, inf)
    pv = chi2(df=dof).sf(lrs)
    return clip(pv, epsilon.super_tiny, 1 - epsilon.tiny)
