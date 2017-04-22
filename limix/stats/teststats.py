import scipy as sp


def empirical_pvalues(xt, x0):
    r"""Function to compute empirical P values.

    Compute empirical P values from the test statistics
    observed on the data and the null test statistics
    (from permutations, parametric bootstraps, etc).

    Args:
        xt (array-like):
            test statistcs observed on data.
        x0 (array-like):
            null test statistcs.
            The minimum P value that can be
            estimated is ``1./float(len(x0))``.

    Returns:
        pv (ndarray):
            estimated empirical P values.

    Example
    -------

        .. doctest::

            >>> from numpy.random import RandomState
            >>> from limix.stats import empirical_pvalues
            >>> from numpy import set_printoptions
            >>> set_printoptions(4)
            >>> random = RandomState(1)
            >>>
            >>> xt = random.chisquare(1, 1000)
            >>> x0 = random.chisquare(1, 10000)
            >>>
            >>> pv = empirical_pvalues(xt, x0)
            >>>
            >>> print(pv.shape)
            (1000,)
            >>>
            >>> print(pv[:4])
            [ 0.5599  1.      0.8389  0.7975]
    """
    idxt = sp.argsort(xt)[::-1]
    idx0 = sp.argsort(x0)[::-1]
    xts = xt[idxt]
    x0s = x0[idx0]
    it = 0
    i0 = 0
    _count = 0
    count = sp.zeros(xt.shape[0])
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
