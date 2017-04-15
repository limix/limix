def boxcox(X):
    """Gaussianize X using the Box-Cox transformation.

    Each phentoype is brought to a positive schale by first subtracting the
    minimum value and adding 1.
    Then each phenotype is transformed by the Box-Cox transformation.

    Args:
        X (array_like): samples by phenotypes.

    Returns:
        array_like: Box-Cox power transformed array.
    """
    from limix_legacy.utils.preprocess import boxcox as _boxcox

    return _boxcox(X)[0]
