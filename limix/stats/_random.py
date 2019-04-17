def multivariate_normal(random, mean, cov):
    from numpy.linalg import cholesky

    L = cholesky(cov)
    return L @ random.randn(L.shape[0]) + mean
