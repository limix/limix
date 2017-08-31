from __future__ import division

from numpy import diag
from numpy_sugar.linalg import economic_qs

from glimix_core.glmm import GLMM
from glimix_core.lmm import LMM

from .model import QTLModel
from .util import (
    assure_named, covariates_process, kinship_process, phenotype_process,
    print_analysis
)


def scan(G, y, lik, K=None, M=None, verbose=True):
    r"""Single-variant association testing via generalised linear mixed models.

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
    :class:`limix.qtl.model.QTLModel`
        QTL representation.

    Examples
    --------
    .. doctest::

        >>> from numpy import dot, exp, sqrt, ones
        >>> from numpy.random import RandomState
        >>> from pandas import DataFrame, set_option
        >>> from limix.qtl import scan
        >>>
        >>> random = RandomState(1)
        >>> set_option('precision', 5)
        >>>
        >>> n = 30
        >>> p = 3
        >>>
        >>> M = DataFrame(dict(offset=ones(n), age=random.randint(10, 60, n)))
        >>>
        >>> X = random.randn(n, 100)
        >>> K = dot(X, X.T)
        >>>
        >>> candidates = random.randn(n, p)
        >>> candidates = DataFrame(candidates, columns=['rs0', 'rs1', 'rs2'])
        >>>
        >>> y = random.poisson(exp(random.randn(n)))
        >>>
        >>> model = scan(candidates, y, 'poisson', K, M=M, verbose=False)
        >>>
        >>> print(model.variant_pvalues.to_string())
        rs0    0.68333
        rs1    0.28812
        rs2    0.51494
        >>> print(model.variant_effsizes.to_string())
        rs0   -0.08455
        rs1   -0.26729
        rs2   -0.15348
        >>> print(model.variant_effsizes_se.to_string())
        rs0    0.20726
        rs1    0.25162
        rs2    0.23570
        >>> print(model)
        Variants
               effsizes  effsizes_se  pvalues
        count   3.00000      3.00000  3.00000
        mean   -0.16844      0.23153  0.49546
        std     0.09228      0.02247  0.19832
        min    -0.26729      0.20726  0.28812
        25%    -0.21038      0.22148  0.40153
        50%    -0.15348      0.23570  0.51494
        75%    -0.11902      0.24366  0.59913
        max    -0.08455      0.25162  0.68333
        <BLANKLINE>
        Covariate effect sizes for the null model
             age  offset
        -0.00557  0.3953
    """
    lik = lik.lower()

    if verbose:
        print_analysis(lik, "Quantitative trait locus analysis")

    y = phenotype_process(lik, y)

    nsamples = len(G)
    G = assure_named(G)

    mixed = K is not None

    M = covariates_process(M, nsamples)

    K, QS = kinship_process(K, nsamples, verbose)

    if lik == 'normal':
        model = _perform_lmm(y, M, QS, G, mixed, verbose)
    else:
        model = _perform_glmm(y, lik, M, K, QS, G, mixed, verbose)

    if verbose:
        print(model)

    return model


def _perform_lmm(y, M, QS, G, mixed, verbose):
    from pandas import Series

    lmm = LMM(y, M.values, QS)
    if not mixed:
        lmm.delta = 1
        lmm.fix('delta')

    lmm.learn(verbose=verbose)

    null_lml = lmm.lml()

    beta = lmm.beta

    keys = list(M.keys())
    ncov_effsizes = Series(beta, keys)

    flmm = lmm.get_fast_scanner()
    alt_lmls, effsizes = flmm.fast_scan(G.values, verbose=verbose)

    alt_lmls = Series(alt_lmls, list(G.columns))
    effsizes = Series(effsizes, list(G.columns))

    return QTLModel(null_lml, alt_lmls, effsizes, ncov_effsizes)


def _perform_glmm(y, lik, M, K, QS, G, mixed, verbose):
    from pandas import Series
    glmm = GLMM(y, lik, M.values, QS)
    if not mixed:
        glmm.delta = 1
        glmm.fix('delta')
    glmm.feed().maximize(verbose=verbose)

    # extract stuff from glmm
    eta = glmm._site.eta
    tau = glmm._site.tau
    scale = float(glmm.scale)
    delta = float(glmm.delta)

    beta = glmm.beta

    keys = list(M.keys())
    ncov_effsizes = Series(beta, keys)

    # define useful quantities
    mu = eta / tau
    var = 1. / tau
    s2_g = scale * (1 - delta)
    tR = diag(var - var.min() + 1e-4)
    tR += s2_g * K

    lmm = LMM(mu, X=M.values, QS=economic_qs(tR))
    lmm.learn(verbose=verbose)
    null_lml = lmm.lml()
    flmm = lmm.get_fast_scanner()

    alt_lmls, effsizes = flmm.fast_scan(G.values, verbose=verbose)

    alt_lmls = Series(alt_lmls, list(G.keys()))
    effsizes = Series(effsizes, list(G.keys()))

    return QTLModel(null_lml, alt_lmls, effsizes, ncov_effsizes)
