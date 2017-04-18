from numpy import copyto

def gower_norm(K, out=None):
    r"""Perform Gower rescaling of covariance matrix K.

    The rescaled covariance matrix has sample variance of 1.

    Example
    -------

        .. doctest::

            >>> from numpy.random import RandomState
            >>> from limix.stats import gower_norm 
            >>> import scipy as sp
            >>>
            >>> X = RandomState(1).randn(4, 4)
            >>> K = sp.dot(X,X.T)
            >>> Z = sp.random.multivariate_normal(sp.zeros(4), K, 50)
            >>> print sp.mean(Z.var(1))
            2.00178627286
            >>> Kn = gower_norm(K)
            >>> Zn = sp.random.multivariate_normal(sp.zeros(4), Kn, 50)
            >>> print sp.mean(Zn.var(1))
            0.748747252259
    """


    c = (K.shape[0] - 1) / (K.trace() - K.mean(0).sum())
    if out is None:
        return c * K

    copyto(out, K)
    out *= c
