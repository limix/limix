from __future__ import division

from numpy import ascontiguousarray, copy, ones, var
from numpy_sugar.linalg import economic_qs

from glimix_core.glmm import GLMMExpFam


def estimate(pheno, lik, K, covs=None, verbose=True):
    r"""Estimate the so-called narrow-sense heritability.

    It supports Normal, Bernoulli, Binomial, and Poisson phenotypes.
    Let :math:`N` be the sample size and :math:`S` the number of covariates.

    Parameters
    ----------
    pheno : tuple, array_like
        Phenotype. Dimensions :math:`N\\times 0`.
    lik : {'normal', 'bernoulli', 'binomial', 'poisson'}
        Likelihood name.
    K : array_like
        Kinship matrix. Dimensions :math:`N\\times N`.
    covs : array_like
        Covariates. Default is an offset. Dimensions :math:`N\\times S`.

    Returns
    -------
    float
        Estimated heritability.

    Examples
    --------
    .. doctest::

        >>> from numpy import dot, exp, sqrt
        >>> from numpy.random import RandomState
        >>> from limix.heritability import estimate
        >>>
        >>> random = RandomState(0)
        >>>
        >>> G = random.randn(50, 100)
        >>> K = dot(G, G.T)
        >>> z = dot(G, random.randn(100)) / sqrt(100)
        >>> y = random.poisson(exp(z))
        >>>
        >>> print('%.2f' % estimate(y, 'poisson', K, verbose=False))
        0.70
    """

    K = _background_standardize(K)
    QS = economic_qs(K)

    lik = lik.lower()

    if lik == "binomial":
        p = len(pheno[0])
    else:
        p = len(pheno)

    if covs is None:
        covs = ones((p, 1))

    glmm = GLMMExpFam(pheno, lik, covs, QS)
    glmm.feed().maximize(verbose=verbose)

    g = glmm.scale * (1 - glmm.delta)
    e = glmm.scale * glmm.delta
    h2 = g / (var(glmm.mean()) + g + e)

    return h2


def _background_standardize(K):
    from ..stats.kinship import gower_norm

    K = copy(K, "C")
    K = ascontiguousarray(K, dtype=float)
    gower_norm(K, K)
    K /= K.diagonal()

    return K
