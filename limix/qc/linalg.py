def remove_dependent_cols(X, tol=1e-6, verbose=False):
    r"""Remove dependent columns.

    Return a matrix with dependent columns removed.

    Parameters
    ----------
    X : array_like
        Matrix to might have dependent columns.

    Returns
    -------
    array_like
        Full column rank matrix.
    """
    from scipy.linalg import qr
    from numpy import abs as npabs
    from numpy import any as npany
    from numpy import where

    R = qr(X, mode="r")[0][: X.shape[1], :]
    I = npabs(R.diagonal()) > tol
    if npany(~I) and verbose:
        msg = "Columns " + str(where(~I)[0])
        print(msg + " have been removed because linear dependence")
        R = X[:, I]
    else:
        R = X.copy()
    return R
