def regress_out(Y, X, return_b=False):
    """
    regresses out X from Y
    """
    Xd = la.pinv(X)
    b = Xd.dot(Y)
    Y_out = Y - X.dot(b)
    if return_b:
        return Y_out, b
    else:
        return Y_out
