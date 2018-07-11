def pca(X, ncomp):
    r"""Principal component analysis.

    Args:
        X (array_like): data.
        ncomp (int): number of components.

    Returns:
        dict: dict containing:
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
            >>> from limix.stats import pca
            >>> from numpy import set_printoptions
            >>> set_printoptions(4)
            >>>
            >>> X = RandomState(1).randn(4, 5)
            >>> result = pca(X, ncomp=2)
            >>> print(result['components']) # doctest: +SKIP
            [[-0.7502  0.5835 -0.0797  0.1957 -0.2285]
             [ 0.4884  0.7227  0.0197 -0.4616 -0.1603]]
            >>> print(result['explained_variance'])
            [6.4466 0.5145]
            >>> print(result['explained_variance_ratio'])
            [0.9205 0.0735]
    """
    from sklearn.decomposition import PCA

    pca = PCA(n_components=ncomp)
    pca.fit(X)

    return dict(
        components=pca.components_,
        explained_variance=pca.explained_variance_,
        explained_variance_ratio=pca.explained_variance_ratio_,
    )
