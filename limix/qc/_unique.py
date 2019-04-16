def unique_variants(X):
    """
    Filters out variants with the same genetic profile.

    Parameters
    ----------
    X : array_like
        Samples-by-variants matrix of genotype values.

    Returns
    -------
    genotype : ndarray
        Genotype array with unique variants.

    Example
    -------
    .. doctest::

        >>> from numpy.random import RandomState
        >>> from numpy import kron, ones, sort
        >>> from limix.qc import unique_variants
        >>> random = RandomState(1)
        >>>
        >>> N = 4
        >>> X = kron(random.randn(N, 3) < 0, ones((1, 2)))
        >>>
        >>> print(X)
        [[0. 0. 1. 1. 1. 1.]
         [1. 1. 0. 0. 1. 1.]
         [0. 0. 1. 1. 0. 0.]
         [1. 1. 0. 0. 1. 1.]]
        >>>
        >>> print(unique_variants(X))
        [[0. 1. 1.]
         [1. 1. 0.]
         [0. 0. 1.]
         [1. 1. 0.]]
    """
    from limix._bits import dask
    import numpy as np

    if dask.is_array(X):
        import dask.array as da

        dot = da.dot
        unique = da.unique
    else:

        dot = np.dot
        unique = np.unique

    u = np.random.rand(X.shape[0])
    i = unique(dot(u, X), return_index=True)[1]
    return X[:, i]
