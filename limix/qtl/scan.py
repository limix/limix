from __future__ import division

from numpy import all as npall
from numpy import diag, eye, isfinite, ones
from numpy_sugar.linalg import economic_qs

from glimix_core.glmm import GLMM, GLMMExpFam
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
        >>> from pandas import DataFrame, option_context
        >>> from limix.qtl import scan
        >>>
        >>> random = RandomState(1)
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
        rs0    0.683330
        rs1    0.288122
        rs2    0.514938
        >>> print(model.variant_effsizes.to_string())
        rs0   -0.084546
        rs1   -0.267287
        rs2   -0.153484
        >>> print(model.variant_effsizes_se.to_string())
        rs0    0.207261
        rs1    0.251623
        rs2    0.235705
        >>> print(model)
        Variants
               effsizes  effsizes_se   pvalues
        count  3.000000     3.000000  3.000000
        mean  -0.168439     0.231530  0.495463
        std    0.092284     0.022474  0.198323
        min   -0.267287     0.207261  0.288122
        25%   -0.210385     0.221483  0.401530
        50%   -0.153484     0.235705  0.514938
        75%   -0.119015     0.243664  0.599134
        max   -0.084546     0.251623  0.683330
        <BLANKLINE>
        Covariate effect sizes for the null model
                age    offset
        -0.00556721  0.395267
    """
    lik = lik.lower()

    if verbose:
        print_analysis(lik, "Quantitative trait locus analysis")

    y = phenotype_process(lik, y)

    nsamples = len(G)
    G = assure_named(G)
    if not npall(isfinite(G)):
        raise ValueError("Variant values must be finite.")

    mixed = K is not None

    M = covariates_process(M, nsamples)

    K, QS = kinship_process(K, nsamples, verbose)

    if lik == 'normal':
        model = _perform_lmm(y, M, QS, G, mixed, verbose)
    else:
        model = _perform_glmm_version1(y, lik, M, K, QS, G, mixed, verbose)

    if verbose:
        print(model)

    return model


def _perform_lmm(y, M, QS, G, mixed, verbose):
    from pandas import Series

    lmm = LMM(y, M.values, QS)
    if not mixed:
        lmm.delta = 1
        lmm.fix('delta')

    lmm.fit(verbose=verbose)

    null_lml = lmm.lml()

    beta = lmm.beta

    keys = list(M.keys())
    ncov_effsizes = Series(beta, keys)

    flmm = lmm.get_fast_scanner()
    alt_lmls, effsizes = flmm.fast_scan(G.values, verbose=verbose)

    alt_lmls = Series(alt_lmls, list(G.columns))
    effsizes = Series(effsizes, list(G.columns))

    return QTLModel(null_lml, alt_lmls, effsizes, ncov_effsizes)


def _perform_glmm_version1(y, lik, M, K, QS, G, mixed, verbose):
    from pandas import Series
    glmm = GLMMExpFam(y, lik, M.values, QS)
    if not mixed:
        glmm.delta = 1
        glmm.fix('delta')
    glmm.fit(verbose=verbose)

    eta = glmm.site.eta
    tau = glmm.site.tau
    scale = float(glmm.scale)
    delta = float(glmm.delta)

    beta = glmm.beta

    keys = list(M.keys())
    ncov_effsizes = Series(beta, keys)

    mu = eta / tau
    var = 1. / tau
    s2_g = scale * (1 - delta)
    tR = diag(var - var.min() + 1e-4)
    tR += s2_g * K

    lmm = LMM(mu, X=M.values, QS=economic_qs(tR))
    lmm.fit(verbose=verbose)
    null_lml = lmm.lml()
    flmm = lmm.get_fast_scanner()

    alt_lmls, effsizes = flmm.fast_scan(G.values, verbose=verbose)

    alt_lmls = Series(alt_lmls, list(G.keys()))
    effsizes = Series(effsizes, list(G.keys()))

    return QTLModel(null_lml, alt_lmls, effsizes, ncov_effsizes)


def _perform_glmm_version2(y, lik, M, K, QS, G, mixed, verbose):
    pass


def _perform_glmm_version3(y, lik, M, K, QS, G, mixed, verbose):
    pass


def _perform_glmm_version4(y, lik, M, K, QS, G, mixed, verbose):
    pass


def _kinship_process(K, nsamples, verbose):

    if K is None:
        K = eye(nsamples)

    K = asarray(K)

    ah = array_hash(K)
    nvalid = _cache['K']['K'] is None or (_cache['K']['hash'] != ah)

    if nvalid:

        if not all(isfinite(K)):
            msg = "One or more values of the provided covariance matrix "
            msg += "is not finite."
            raise ValueError(msg)

        K = gower_norm(K)

        desc = "Eigen decomposition of the covariance matrix..."
        with Timer(desc=desc, disable=not verbose):
            QS = economic_qs(K)
            _cache['K']['hash'] = ah
            _cache['K']['QS'] = QS
            _cache['K']['K'] = K
    else:
        QS = _cache['K']['QS']
        K = _cache['K']['K']

    return K, QS


def _phenotype_process(y):
    if isinstance(y, (tuple, list)):
        y = tuple([asarray(p, float) for p in y])
    else:
        y = asarray(y, float)

    if not all(isfinite(y)):
        msg = "One or more values of the provided phenotype "
        msg += "is not finite."
        raise ValueError(msg)
    return y


def _covariates_process(M, nsamples):
    from pandas import DataFrame

    if M is None:
        M = DataFrame({'offset': ones(nsamples, float)})
    else:
        M = assure_named(M, nsamples)

    if not all(isfinite(M.values)):
        msg = "One or more values of the provided covariates "
        msg += "is not finite."
        raise ValueError(msg)

    return M
