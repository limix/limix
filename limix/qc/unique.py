def unique_variants(X):
    r"""Filters out variants with the same genetic profile.

    Parameters
    ----------
    X : ndarray
        (`N`, `S`) ndarray of genotype values for `N` individuals and `S`
        variants.

    Returns
    -------
    ndarray
        Genotype array with unique variants.

    Examples
    --------
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
        >>> idx = unique_variants(X)
        >>>
        >>> print(X[:, sort(idx)])
        [[0. 1. 1.]
         [1. 0. 1.]
         [0. 1. 0.]
         [1. 0. 1.]]
    """
    import dask.array as da
    import numpy as np

    if isinstance(X, da.Array):
        dot = da.dot
        unique = _dask_unique
    else:
        dot = np.dot
        unique = np.unique

    u = np.random.rand(X.shape[0])
    return unique(dot(u, X), return_index=True)[1]


def _dask_unique(x, return_index=True):
    from dask.array.core import Array
    from dask import sharedict
    from numpy import cumsum, concatenate, unique
    from numpy.testing import assert_

    assert_(return_index)

    name = "unique-" + x.name

    def _unique(x):

        return unique(x, return_index=return_index)

    dsk = dict(((name, i), (_unique, key)) for i, key in enumerate(x._keys()))
    parts = Array._get(sharedict.merge((name, dsk), x.dask), list(dsk.keys()))

    arrs = [a[0] for a in parts]

    chunks = x.chunks[0]
    offset = cumsum((0,) + chunks)[:-1]

    idxs = [parts[i][1] + offset[i] for i in range(len(parts))]

    arr = concatenate(arrs)
    idx = concatenate(idxs)

    u, i = unique(arr, return_index=True)
    return u, idx[i]
