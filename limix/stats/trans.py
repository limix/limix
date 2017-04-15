def boxcox(X):
    """Gaussianize X using the Box-Cox transformation.

    Each phentoype is brought to a positive schale by first subtracting the
    minimum value and adding 1.
    Then each phenotype is transformed by the Box-Cox transformation.

    Args:
        X (array_like): samples by phenotypes.

    Returns:
        array_like: Box-Cox power transformed array.

    Example
    -------

        .. doctest::

            >>> from numpy.random import RandomState
            >>> from limix.stats import boxcox
            >>>
            >>> random = RandomState(0)
            >>> X = random.randn(5, 2)
            >>>
            >>> print(boxcox(X))
            [[ 2.71356378  0.95441669]
             [ 1.38440507  1.69459001]
             [ 2.90661679  0.        ]
             [ 1.34068423  0.64396462]
             [ 0.          0.9597255 ]]
    """
    from limix_legacy.utils.preprocess import boxcox as _boxcox

    return _boxcox(X)[0]
