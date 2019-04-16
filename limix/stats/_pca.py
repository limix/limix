def pca(X, ncomp):
    """
    Principal component analysis.

    Parameters
    ----------
    X : array_like
        Samples-by-dimensions array.
    ncomp : int
        Number of components.

    Returns
    -------
    components : ndarray
        First components ordered by explained variance.
    explained_variance : ndarray
        Explained variance.
    explained_variance_ratio : ndarray
        Percentage of variance explained.

    Examples
    --------
    .. doctest::

        >>> from numpy import round
        >>> from numpy.random import RandomState
        >>> from limix.stats import pca
        >>>
        >>> X = RandomState(1).randn(4, 5)
        >>> r = pca(X, ncomp=2)
        >>> r['components']
        array([[-0.75015369,  0.58346541, -0.07973564,  0.19565682, -0.22846925],
               [ 0.48842769,  0.72267548,  0.01968344, -0.46161623, -0.16031708]])
        >>> r['explained_variance'] # doctest: +FLOAT_CMP
        array([6.44655993, 0.51454938])
        >>> r['explained_variance_ratio'] # doctest: +FLOAT_CMP
        array([0.92049553, 0.07347181])
    """
    from sklearn.decomposition import PCA
    from numpy import asarray

    X = asarray(X, float)
    pca = PCA(n_components=ncomp)
    pca.fit(X)

    return dict(
        components=pca.components_,
        explained_variance=pca.explained_variance_,
        explained_variance_ratio=pca.explained_variance_ratio_,
    )
