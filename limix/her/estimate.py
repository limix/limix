from __future__ import division

from numpy import pi, var

from glimix_core.glmm import GLMMExpFam
from glimix_core.lmm import LMM
from numpy_sugar.linalg import economic_qs

from ..dataframe import normalise_dataset
from ..fprint import oprint
from ..likelihood import normalise_extreme_values
from ..qc import normalise_covariance


def estimate(y, lik, K, M=None, verbose=True):
    r"""Estimate the so-called narrow-sense heritability.

    It supports Normal, Bernoulli, Binomial, and Poisson phenotypes.
    Let :math:`N` be the sample size and :math:`S` the number of covariates.

    Parameters
    ----------
    y : array_like
        Either a tuple of two arrays of `N` individuals each (Binomial
        phenotypes) or an array of `N` individuals (Normal, Poisson, or
        Bernoulli phenotypes). It does not support missing values yet.
    lik : 'normal', 'bernoulli', 'binomial', 'poisson'
        Sample likelihood describing the residual distribution.
    K : array_like
        `N` by `N` covariance matrix (e.g., kinship coefficients).
    M : array_like, optional
        `N` individuals by `D` covariates.
        By default, ``M`` is a (`N`, `1`) array of ones.
    verbose : bool, optional
        if ``True``, details such as runtime are displayed.

    Returns
    -------
    float
        Estimated heritability.

    Examples
    --------
    .. doctest::

        >>> from numpy import dot, exp, sqrt
        >>> from numpy.random import RandomState
        >>> from limix.her import estimate
        >>>
        >>> random = RandomState(0)
        >>>
        >>> G = random.randn(150, 200) / sqrt(200)
        >>> K = dot(G, G.T)
        >>> z = dot(G, random.randn(200)) + random.randn(150)
        >>> y = random.poisson(exp(z))
        >>>
        >>> print('%.2f' % estimate(y, 'poisson', K, verbose=False))
        0.18
    """
    if verbose:
        lik_name = lik.lower()
        lik_name = lik_name[0].upper() + lik_name[1:]
        analysis_name = "Heritability estimation"
        oprint("*** %s using %s-GLMM ***" % (analysis_name, lik_name))

    if K is not None:
        K = normalise_covariance(K)

    data = normalise_dataset(y, lik, M=M, K=K)
    y = data['y']
    M = data['M']
    K = data['K']

    y = normalise_extreme_values(y, lik)

    if K is not None:
        QS = economic_qs(K)
    else:
        QS = None

    if lik == 'normal':
        method = LMM(y.values, M.values, QS)
        method.fit(verbose=verbose)
    else:
        method = GLMMExpFam(y.values, lik, M.values, QS, n_int=500)
        method.fit(verbose=verbose, factr=1e6, pgtol=1e-3)

    g = method.scale * (1 - method.delta)
    e = method.scale * method.delta
    if lik == 'bernoulli':
        e += pi * pi / 3

    if lik == 'normal':
        v = method.fixed_effects_variance
    else:
        v = var(method.mean())

    return g / (v + g + e)
