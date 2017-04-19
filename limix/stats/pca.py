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

    Example
    -------

        .. doctest::

            >>> from numpy.random import RandomState
            >>> from limix.stats import pca
            >>>
            >>> X = RandomState(1).randn(4, 5)
            >>> result = pca(X, ncomp=2)
            >>> print(result['components'])
            [[-0.75015369  0.58346541 -0.07973564  0.19565682 -0.22846925]
             [ 0.48842769  0.72267548  0.01968344 -0.46161623 -0.16031708]]
            >>> print(result['explained_variance'])
            [ 4.83491994  0.38591204]
            >>> print(result['explained_variance_ratio'])
            [ 0.92049553  0.07347181]
    """
    from sklearn.decomposition import PCA

    pca = PCA(n_components=ncomp)
    pca.fit(X)

    return dict(
        components=pca.components_,
        explained_variance=pca.explained_variance_,
        explained_variance_ratio=pca.explained_variance_ratio_)
