from numpy import all as npall
from numpy import concatenate, isfinite, newaxis
from tqdm import tqdm

from glimix_core.glmm import GLMM
from glimix_core.lmm import LMM

from .model import IQTLModel
from .util import (
    assure_named, covariates_process, kinship_process, phenotype_process,
    print_analysis
)


def iscan(G, y, lik, inter, K=None, M=None, verbose=True):
    r"""Interaction single-variant association testing via mixed models.

    It supports Normal (linear mixed model), Bernoulli, Binomial, and Poisson
    residual errors, defined by ``lik``.
    The columns of ``G`` define the candidates to be tested for association
    with the phenotype ``y``.
    The covariance matrix is set by ``K``.
    If not provided, or set to ``None``, the generalised linear model
    without random effects is assumed.
    The covariates can be set via the parameter ``M``.
    We recommend to always provide a column of ones in the case

    Parameters
    ----------
    G : array_like
        `n` individuals by `s` candidate markers.
    y : tuple, array_like
        Either a tuple of two arrays of `n` individuals each (Binomial
        phenotypes) or an array of `n` individuals (Normal, Poisson, or
        Bernoulli phenotypes). It does not support missing values yet.
    lik : {'normal', 'bernoulli', 'binomial', 'poisson'}
        Sample likelihood describing the residual distribution.
    inter : array_like
        `n` individuals by `i` interaction factors.
    K : array_like, optional
        `n` by `n` covariance matrix (e.g., kinship coefficients).
        Set to ``None`` for a generalised linear model without random effects.
        Defaults to ``None``.
    M : array_like, optional
        `n` individuals by `d` covariates.
        By default, ``M`` is an `n` by `1` matrix of ones.
    verbose : bool, optional
        if ``True``, details such as runtime are displayed.

    Returns
    -------
    :class:`limix.qtl.model.IQTLModel`
        Interaction QTL representation.

    Examples
    --------
    .. doctest::

        >>> from numpy import dot, zeros, asarray
        >>> from numpy.random import RandomState
        >>> from numpy.testing import assert_allclose
        >>>
        >>> from limix.qtl import iscan
        >>> from pandas import DataFrame, option_context
        >>>
        >>> random = RandomState(0)
        >>> nsamples = 50
        >>>
        >>> X = random.randn(nsamples, 10)
        >>> G = random.randn(nsamples, 100)
        >>> K = dot(G, G.T)
        >>> ntrials = random.randint(1, 100, nsamples)
        >>> z = dot(G, random.randn(100)) / 10
        >>>
        >>> successes = zeros(len(ntrials), int)
        >>> for i in range(len(ntrials)):
        ...     for j in range(ntrials[i]):
        ...         successes[i] += int(z[i] + 0.5 * random.randn() > 0)
        >>>
        >>> y = successes / asarray(ntrials, float)
        >>>
        >>> inter = random.randn(nsamples, 3)
        >>>
        >>> index = ['sample%02d' % i for i in range(X.shape[0])]
        >>> cols = ['SNP%02d' % i for i in range(X.shape[1])]
        >>> X = DataFrame(data=X, index=index, columns=cols)
        >>>
        >>> cols = ['inter%02d' % i for i in range(inter.shape[1])]
        >>> inter = DataFrame(data=inter, index=index, columns=cols)
        >>>
        >>> model = iscan(X, y, 'normal', inter, K, verbose=False)
        >>>
        >>> with option_context('precision', 5):
        ...     print(model.variant_pvalues)
               inter00  inter01  inter02
        SNP00  0.81180  0.63035  0.61240
        SNP01  0.02847  0.64437  0.82671
        SNP02  0.56817  0.72882  0.23928
        SNP03  0.53793  0.64628  0.86144
        SNP04  0.13858  0.39475  0.28650
        SNP05  0.06722  0.56295  0.39859
        SNP06  0.12739  0.62219  0.68084
        SNP07  0.32834  0.96894  0.67628
        SNP08  0.28341  0.29361  0.56248
        SNP09  0.64945  0.67185  0.76600
    """
    lik = lik.lower()

    if verbose:
        print_analysis(lik, "Interaction QTL analysis")

    y = phenotype_process(lik, y)

    nsamples = len(G)
    G = assure_named(G)
    if not npall(isfinite(G)):
        raise ValueError("Variant values must be finite.")

    inter = assure_named(inter)
    if not npall(isfinite(inter)):
        raise ValueError("Interaction values must be finite.")

    mixed = K is not None

    M = covariates_process(M, nsamples)

    K, QS = kinship_process(K, nsamples, verbose)

    if lik == 'normal':
        model = _perform_lmm(y, M, QS, G, inter, mixed, verbose)
    else:
        raise NotImplementedError

    if verbose:
        print(model)

    return model


def _perform_lmm(y, M, QS, G, inter, mixed, verbose):
    from pandas import DataFrame, Series

    alt_lmls = dict()
    effsizes = dict()
    ncov_effsizes = dict()
    null_lmls = []
    interv = inter.values

    for c in tqdm(G.columns, disable=not verbose):
        g = G[c].values[:, newaxis]
        X1 = g * interv

        covariates = concatenate((M.values, g), axis=1)
        lmm = LMM(y, covariates, QS)
        if not mixed:
            lmm.delta = 1
            lmm.fix('delta')

        lmm.fit(verbose=False)

        null_lmls.append(lmm.lml())

        ncov_effsizes[c] = lmm.beta

        flmm = lmm.get_fast_scanner()
        alt_lmls[c], effsizes[c] = flmm.fast_scan(X1, verbose=False)

    alt_lmls = DataFrame(data=alt_lmls, index=inter.columns).transpose()
    effsizes = DataFrame(data=effsizes, index=inter.columns).transpose()

    index = list(M.columns) + ["variant"]
    ncov_effsizes = DataFrame(data=ncov_effsizes, index=index).transpose()
    null_lml = Series(null_lmls, G.columns)

    return IQTLModel(null_lml, alt_lmls, effsizes, ncov_effsizes)
