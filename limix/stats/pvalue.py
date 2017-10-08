from __future__ import division


# TODO: add test
def multipletests(pvals, alpha=0.05, method='hs', is_sorted=False):
    r'''Test results and p-value correction for multiple tests.

    Parameters
    ----------
    pvals : array_like
        Uncorrected p-values.
    alpha : float
        FWER, family-wise error rate, e.g. ``0.1``.
    method : string
        Method used for testing and adjustment of pvalues. Can be either the
        full name or initial letters. Available methods are ::

        `bonferroni` : one-step correction
        `sidak` : one-step correction
        `holm-sidak` : step down method using Sidak adjustments
        `holm` : step-down method using Bonferroni adjustments
        `simes-hochberg` : step-up method  (independent)
        `hommel` : closed method based on Simes tests (non-negative)
        `fdr_bh` : Benjamini/Hochberg  (non-negative)
        `fdr_by` : Benjamini/Yekutieli (negative)
        `fdr_tsbh` : two stage fdr correction (non-negative)
        `fdr_tsbky` : two stage fdr correction (non-negative)

    is_sorted : bool
        If ``False`` (default), the p_values will be sorted, but the corrected
        pvalues are in the original order. If ``True``, then it assumed that
        the pvalues are already sorted in ascending order.

    Returns
    -------
    reject : array_like, boolean
        ``true`` for hypothesis that can be rejected for given alpha.
    pvals_corrected : array_like
        p-values corrected for multiple tests.
    alphacSidak: float
        corrected alpha for Sidak method.
    alphacBonf: float
        corrected alpha for Bonferroni method.

    Notes
    -----
    This is a wrapper to a function from the `statsmodels`_ package.

    .. _statsmodels: http://www.statsmodels.org
    '''

    from statsmodels.sandbox.stats.multicomp import multipletests as mt
    return mt(pvals, alpha=alpha, method=method, is_sorted=is_sorted)
