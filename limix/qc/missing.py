from __future__ import division


def count_missingness(X):
    r"""Count the number of missing values per column.

    Returns
    -------
    array_like
        Number of missing values per column.
    """
    import dask.array as da
    from numpy import isnan

    if isinstance(X, da.Array):
        return da.isnan(X).sum(axis=0).compute()

    return isnan(X).sum(axis=0)
