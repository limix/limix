def boxcox(X):
    r"""Gaussianize X using the Box-Cox transformation.

    Each phentoype is brought to a positive schale by first subtracting the
    minimum value and adding 1.
    Then each phenotype is transformed by the Box-Cox transformation.

    Args:
        X (array_like): samples by phenotypes.

    Returns:
        array_like: Box-Cox power transformed array.

    Examples
    --------

        .. doctest::

            >>> from numpy.random import RandomState
            >>> from limix.stats import boxcox
            >>> from numpy import set_printoptions
            >>> set_printoptions(4)
            >>>
            >>> random = RandomState(0)
            >>> X = random.randn(5, 2)
            >>>
            >>> print(boxcox(X))
            [[2.7136 0.9544]
             [1.3844 1.6946]
             [2.9066 0.    ]
             [1.3407 0.644 ]
             [0.     0.9597]]
    """
    from limix.util.preprocess import boxcox as _boxcox

    return _boxcox(X)[0]
