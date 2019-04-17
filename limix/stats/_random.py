def multivariate_normal(random, mean, cov):
    """
    Draw random samples from a multivariate normal distribution.

    Parameters
    ----------
    random : np.random.RandomState instance
        Random state.
    mean : array_like
        Mean of the n-dimensional distribution.
    cov : array_like
        Covariance matrix of the distribution. It must be symmetric and
        positive-definite for proper sampling.

    Returns
    -------
    out : ndarray
        The drawn sample.
    """
    from numpy.linalg import cholesky

    L = cholesky(cov)
    return L @ random.randn(L.shape[0]) + mean
