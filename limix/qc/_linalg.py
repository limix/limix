def remove_dependent_cols(X, tol=1.49e-08):
    """
    Remove dependent columns.

    Return a matrix with dependent columns removed.

    Parameters
    ----------
    X : array_like
        Matrix to might have dependent columns.
    tol : float, optional
        Threshold above which columns are considered dependents.

    Returns
    -------
    rank : ndarray
        Full column rank matrix.
    """
    from scipy.linalg import qr
    from numpy import abs, asarray
    from numpy import concatenate, full

    X = asarray(X)
    if X.shape[1] == 0:
        return X

    r = max(X.shape[1] - X.shape[0], 0)
    i = abs(qr(X, mode="r")[0].diagonal()) > tol
    i = concatenate((i, full(r, True, bool)))
    x = X[:, i]
    if x.shape[1] == X.shape[1]:
        return X
    return remove_dependent_cols(x, tol)
