def multipletests(pvals, alpha=0.05, method="hs", is_sorted=False):
    """
    Test results and p-value correction for multiple tests.

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
    reject : ndarray, boolean
        ``True`` for hypothesis that can be rejected for given alpha.
    pvals_corrected : ndarray
        P-values corrected for multiple tests.
    alphacSidak : float
        Corrected alpha for Sidak method.
    alphacBonf : float
        Corrected alpha for Bonferroni method.

    Notes
    -----
    This is a wrapper around a function from the `statsmodels`_ package.

    .. _statsmodels: http://www.statsmodels.org
    """

    from statsmodels.sandbox.stats.multicomp import multipletests as mt

    return mt(pvals, alpha=alpha, method=method, is_sorted=is_sorted)


def empirical_pvalues(xt, x0):
    """
    Function to compute empirical p-values.

    Compute empirical p-values from the test statistics
    observed on the data and the null test statistics
    (from permutations, parametric bootstraps, etc).

    Parameters
    ----------
    xt : array_like
        Test statistcs observed on data.
    x0 : array_like
        Null test statistcs. The minimum p-value that can be
        estimated is ``1./len(x0)``.

    Returns
    -------
    pvalues : ndarray
        Estimated empirical p-values.

    Examples
    --------
    .. doctest::

        >>> from numpy.random import RandomState
        >>> from limix.stats import empirical_pvalues
        >>>
        >>> random = RandomState(1)
        >>> x0 = random.chisquare(1, 5)
        >>> x1 = random.chisquare(1, 10000)
        >>>
        >>> empirical_pvalues(x0, x1) # doctest: +FLOAT_CMP
        array([0.563 , 1.    , 0.839 , 0.7982, 0.5803])
    """
    from numpy import argsort, zeros, asarray

    xt = asarray(xt, float)
    x0 = asarray(x0, float)

    idxt = argsort(xt)[::-1]
    idx0 = argsort(x0)[::-1]
    xts = xt[idxt]
    x0s = x0[idx0]
    it = 0
    i0 = 0
    _count = 0
    count = zeros(xt.shape[0])
    while True:
        if x0s[i0] > xts[it]:
            _count += 1
            i0 += 1
            if i0 == x0.shape[0]:
                count[idxt[it:]] = _count
                break
        else:
            count[idxt[it]] = _count
            it += 1
            if it == xt.shape[0]:
                break
    pv = (count + 1) / float(x0.shape[0])
    pv[pv > 1.0] = 1.0
    return pv
