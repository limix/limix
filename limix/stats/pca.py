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

        >>> from numpy.random import RandomState
        >>> from numpy import array_str
        >>> from limix.stats import pca
        >>>
        >>> X = RandomState(1).randn(4, 5)
        >>> result = pca(X, ncomp=2)
        >>> print(array_str(result['components'], precision=4))
        [[-0.7502  0.5835 -0.0797  0.1957 -0.2285]
         [ 0.4884  0.7227  0.0197 -0.4616 -0.1603]]
        >>> print(array_str(result['explained_variance'], precision=4))
        [6.4466 0.5145]
        >>> print(array_str(result['explained_variance_ratio'], precision=4))
        [0.9205 0.0735]
    """
    from sklearn.decomposition import PCA

    pca = PCA(n_components=ncomp)
    pca.fit(X)

    return dict(
        components=pca.components_,
        explained_variance=pca.explained_variance_,
        explained_variance_ratio=pca.explained_variance_ratio_)
