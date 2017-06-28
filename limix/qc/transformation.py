import scipy.linalg as la
from numpy import isnan, zeros_like
from scipy.stats import boxcox as sp_boxcox
from scipy_sugar.stats import quantile_gaussianize


def boxcox(X):
    r"""Gaussianize X using the Box-Cox transformation.

    Each phentoype is brought to a positive schale by first subtracting the
    minimum value and adding 1.
    Then each phenotype is transformed by the Box-Cox transformation.

    Parameters
    ----------
    X : array_like
        Samples by phenotypes.

    Returns
    -------
    array_like
        Box-Cox power transformed array.

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
        [[ 2.7136  0.9544]
         [ 1.3844  1.6946]
         [ 2.9066  0.    ]
         [ 1.3407  0.644 ]
         [ 0.      0.9597]]
    """
    X_transformed = zeros_like(X)
    for i in range(X.shape[1]):
        i_nan = isnan(X[:, i])
        values = X[~i_nan, i]
        X_transformed[i_nan, i] = X[i_nan, i]
        X_transformed[~i_nan, i], _ = sp_boxcox(values - values.min() + 1.0)
    return X_transformed


def mean_standardize(Y, in_place=False):
    r"""Zero-mean and one-deviation normalization.

    Standardize Y in a way that is robust to missing values

    Parameters
    ----------
    Y : array_like
        Array to be normalized.
    in_place : bool
        Whether to operate in-place.

    Returns
    -------
    array_like
        Normalized array.
    """
    if in_place:
        YY = Y
    else:
        YY = Y.copy()
    for i in range(YY.shape[1]):
        Iok = ~SP.isnan(YY[:, i])
        Ym = YY[Iok, i].mean()
        YY[:, i] -= Ym
        Ys = YY[Iok, i].std()
        YY[:, i] /= Ys
    return YY


def remove_dependent_cols(M, tol=1e-6, display=False):
    """Remove dependent columns.

        Return a matrix with dependent columns removed.

    Parameters
    ----------
        M : array_like

    """
    R = la.qr(M, mode='r')[0][:M.shape[1], :]
    I = (abs(R.diagonal()) > tol)
    if sp.any(~I) and display:
        print(('cols ' + str(sp.where(~I)[0]) +
               ' have been removed because linearly dependent on the others'))
        R = M[:, I]
    else:
        R = M.copy()
    return R


def regress_out(Y, X, return_b=False):
    """
    regresses out X from Y
    """
    Xd = la.pinv(X)
    b = Xd.dot(Y)
    Y_out = Y - X.dot(b)
    if return_b:
        return Y_out, b
    else:
        return Y_out
