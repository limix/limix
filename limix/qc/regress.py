def regress_out(Y, X, return_b=False):
    r"""Regresses out X from Y"""
    from scipy.linalg import pinv
    Xd = la.pinv(X)
    b = Xd.dot(Y)
    Y_out = Y - X.dot(b)
    if return_b:
        return Y_out, b
    else:
        return Y_out
