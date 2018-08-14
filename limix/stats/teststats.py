from numpy import argsort, zeros


def empirical_pvalues(xt, x0):
    r"""Function to compute empirical p-values.

    Compute empirical p-values from the test statistics
    observed on the data and the null test statistics
    (from permutations, parametric bootstraps, etc).

    Parameters
    ----------
    xt : array_like
        Test statistcs observed on data.
    x0 : array_like
        Null test statistcs. The minimum p-value that can be
        estimated is ``1./float(len(x0))``.

    Returns
    -------
    array_like
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
        >>> empirical_pvalues(x0, x1) # doctest: +SKIP
        array([0.563 , 1.    , 0.839 , 0.7982, 0.5803])
    """
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
    pv[pv > 1.] = 1.
    return pv
