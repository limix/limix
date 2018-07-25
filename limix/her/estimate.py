from __future__ import division

from numpy import pi, var

from glimix_core.glmm import GLMMExpFam
from glimix_core.lmm import LMM
from numpy_sugar.linalg import economic_qs

from ..dataframe import normalise_dataset
from ..likelihood import normalise_extreme_values
from ..qc import normalise_covariance


def estimate(y, lik, K, M=None, verbose=True):
    r"""Estimate the so-called narrow-sense heritability.

    It supports Normal, Bernoulli, Probit, Binomial, and Poisson phenotypes.
    Let :math:`N` be the sample size and :math:`S` the number of covariates.

    Parameters
    ----------
    y : array_like
        Either a tuple of two arrays of `N` individuals each (Binomial
        phenotypes) or an array of `N` individuals (Normal, Poisson, or
        Bernoulli phenotypes). If a continuous phenotype is provided (i.e., a Normal
        one), make sure they have been normalised in such a way that its values are
        not extremely large; it might cause numerical errors otherwise. For example,
        by using :func:`limix.qc.mean_standardize` or
        :func:`limix.qc.quantile_gaussianize`.
    lik : 'normal', 'bernoulli', 'probit', binomial', 'poisson'
        Sample likelihood describing the residual distribution.
    K : array_like
        `N` by `N` covariance matrix. It might be, for example, the estimated kinship
        relationship between the individuals. The provided matrix will be normalised
        via the function :func:`limix.qc.normalise_covariance`.
    M : array_like, optional
        `N` individuals by `S` covariates.
        It will create a (`N`, `1`) matrix ``M`` of ones representing the offset
        covariate if ``None`` is passed. If an array is passed, it will used as is.
        Defaults to ``None``.
    verbose : bool, optional
        ``True`` to display progress and summary; ``False`` otherwise.

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
        >>> print('%.3f' % estimate(y, 'poisson', K, verbose=False))  # doctest: +FLOAT_CMP
        0.183

    Notes
    -----
    It will raise a ``ValueError`` exception if non-finite values are passed. Please,
    refer to the :func:`limix.qc.mean_impute` function for missing value imputation.
    """
    if verbose:
        lik_name = lik.lower()
        lik_name = lik_name[0].upper() + lik_name[1:]
        analysis_name = "Heritability estimation"
        print("*** {} using {}-GLMM ***".format(analysis_name, lik_name))

    if K is not None:
        K = normalise_covariance(K)

    data = normalise_dataset(y, lik, M=M, K=K)
    y = data["y"]
    M = data["M"]
    K = data["K"]

    y = normalise_extreme_values(y, lik)

    if K is not None:
        QS = economic_qs(K)
    else:
        QS = None

    if lik == "normal":
        method = LMM(y.values, M.values, QS)
        method.fit(verbose=verbose)
    else:
        method = GLMMExpFam(y.values, lik, M.values, QS, n_int=500)
        method.fit(verbose=verbose, factr=1e6, pgtol=1e-3)

    g = method.scale * (1 - method.delta)
    e = method.scale * method.delta
    if lik == "bernoulli":
        e += pi * pi / 3

    if lik == "normal":
        v = method.fixed_effects_variance
    else:
        v = var(method.mean())

    return g / (v + g + e)
