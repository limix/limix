from __future__ import division


def compute_maf(X):
    r"""Compute minor allele frequencies.

    It assumes that ``X`` encodes 0, 1, and 2 representing the number
    of alleles (or dosage), or ``NaN`` to represent missing values.

    Parameters
    ----------
    X : array_like
        Genotype matrix.

    Returns
    -------
    array_like
        Minor allele frequencies.

    Examples
    --------
    .. doctest::

        >>> from numpy.random import RandomState
        >>> from limix.qc import compute_maf
        >>>
        >>> random = RandomState(0)
        >>> X = random.randint(0, 3, size=(100, 10))
        >>>
        >>> print(compute_maf(X)) # doctest: +FLOAT_CMP
        [0.49  0.49  0.445 0.495 0.5   0.45  0.48  0.48  0.47  0.435]
    """
    import dask.array as da
    from numpy import isnan, logical_not, minimum, nansum

    if isinstance(X, da.Array):
        s0 = da.nansum(X, axis=0).compute()
        denom = 2 * (X.shape[0] - da.isnan(X).sum(axis=0)).compute()
    else:
        s0 = nansum(X, axis=0)
        denom = 2 * logical_not(isnan(X)).sum(axis=0)
    s0 = s0 / denom
    s1 = 1 - s0
    return minimum(s0, s1)
