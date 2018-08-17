# TODO: normalise this documentation
def pca(X, ncomp):
    r"""Principal component analysis.

    Parameters
    ----------
    X : array_like
        Data.
    ncomp : int
        Number of components.

    Returns
    -------
    dict
        - **components** (*array_like*):
          first components ordered by explained variance.
        - **explained_variance** (*array_like*):
          explained variance.
        - **explained_variance_ratio** (*array_like*):
          percentage of variance explained.

    Examples
    --------
    .. doctest::

        >>> from numpy import round
        >>> from numpy.random import RandomState
        >>> from limix.stats import pca
        >>>
        >>> X = RandomState(1).randn(4, 5)
        >>> r = pca(X, ncomp=2)
        >>> round(r['components'], 2)
        array([[-0.75,  0.58, -0.08,  0.2 , -0.23],
               [ 0.49,  0.72,  0.02, -0.46, -0.16]])
        >>> round(r['explained_variance'], 4) # doctest: +FLOAT_CMP
        array([ 6.4466,  0.5145])
        >>> round(r['explained_variance_ratio'], 4) # doctest: +FLOAT_CMP
        array([ 0.9205,  0.0735])
    """
    from sklearn.decomposition import PCA

    pca = PCA(n_components=ncomp)
    pca.fit(X)

    return dict(
        components=pca.components_,
        explained_variance=pca.explained_variance_,
        explained_variance_ratio=pca.explained_variance_ratio_,
    )
