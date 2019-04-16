def count_missingness(X):
    """
    Count the number of missing values per column.

    Parameters
    ----------
    X : array_like
        Matrix.

    Returns
    -------
    count : ndarray
        Number of missing values per column.
    """
    import dask.array as da
    from numpy import isnan

    if isinstance(X, da.Array):
        return da.isnan(X).sum(axis=0).compute()

    return isnan(X).sum(axis=0)
