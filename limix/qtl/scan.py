from __future__ import division

from numpy import all as npall
from numpy import asarray, isfinite, ones
from pandas import DataFrame

from glimix_core.glmm import GLMMExpFam, GLMMNormal
from glimix_core.lmm import LMM

from ..nice_arrays import (
    covariates_process, default_covariates_index, kinship_process,
    normalise_candidates_matrix, normalise_covariates_matrix,
    normalise_kinship_matrix, normalise_phenotype_matrix, phenotype_process)
from ..dataframe import infer_samples_index, make_sure_phenotype_dataframe
from ..dataframe import make_sure_covariates_dataframe
from ..dataframe import make_sure_candidates_dataframe
from ..dataframe import make_sure_kinship_dataframe
from ..likelihood import normalise_extreme_values
from .model import QTLModel
from .util import print_analysis


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
        >>> from pandas import DataFrame
        >>> from limix.qtl import scan
        >>>
        >>> random = RandomState(1)
        >>>
        >>> n = 30
        >>> p = 3
        >>> samples_index = range(n)
        >>>
        >>> M = DataFrame(dict(offset=ones(n), age=random.randint(10, 60, n)))
        >>> M.index = samples_index
        >>>
        >>> X = random.randn(n, 100)
        >>> K = dot(X, X.T)
        >>>
        >>> candidates = random.randn(n, p)
        >>> candidates = DataFrame(candidates, index=samples_index,
        ...                                    columns=['rs0', 'rs1', 'rs2'])
        >>>
        >>> y = random.poisson(exp(random.randn(n)))
        >>>
        >>> model = scan(candidates, y, 'poisson', K, M=M, verbose=False)
        >>>
        >>> model.variant_pvalues.round(2)
        rs0    0.55
        rs1    0.22
        rs2    0.55
        dtype: float64
        >>> model.variant_effsizes.round(2)
        rs0   -0.13
        rs1   -0.32
        rs2   -0.14
        dtype: float64
        >>> model.variant_effsizes_se.round(2)
        rs0    0.22
        rs1    0.26
        rs2    0.24
        dtype: float64
        >>> model
        Variants
               effsizes  effsizes_se   pvalues
        count  3.000...     3.000...  3.000...
        mean  -0.196...     0.239...  0.441...
        std    0.102...     0.017...  0.193...
        min   -0.315...     0.221...  0.218...
        25%   -0.229...     0.231...  0.385...
        50%   -0.143...     0.242...  0.552...
        75%   -0.137...     0.249...  0.553...
        max   -0.130...     0.256...  0.554...
        <BLANKLINE>
        Covariate effect sizes for the null model
              age   offset
        -0.005...  0.39...
    """
    lik = lik.lower()

    if verbose:
        print_analysis(lik, "Quantitative trait locus analysis")

    nsamples = len(G)
    if isinstance(y, (tuple, list)):
        y = asarray(y, float).T

    if M is None:
        M = ones((nsamples, 1))

    arrs = [a for a in [y, G, M, K] if a is not None]
    samples_index = infer_samples_index(arrs)
    if len(samples_index) == 0:
        raise ValueError("Could not infer an index for samples."
                         " Please, check if the passed arrays are in order.")

    y = make_sure_phenotype_dataframe(y, samples_index)
    G = make_sure_candidates_dataframe(G, samples_index)
    M = make_sure_covariates_dataframe(M, samples_index)

    if K is not None:
        K = make_sure_kinship_dataframe(K, samples_index)

    y = normalise_extreme_values(y, lik)
    if K is not None:
        K = normalise_kinship_matrix(K)

    # import pdb
    # pdb.set_trace()
    # G = normalise_candidates_matrix(G)
    # # M = covariates_process(M, nsamples)
    # import pdb
    # pdb.set_trace()
    # # M = normalise_covariates_matrix(M)

    # y = phenotype_process(lik, y)

    if not npall(isfinite(G)):
        raise ValueError("Candidate values must be finite.")

    mixed = K is not None

    indices = _intersect_indices(G, y, K, M)

    M = M.loc[indices]
    G = G.loc[indices]
    y = y.loc[indices]

    if K is not None:
        K = K.loc[indices, :].loc[:, indices]

    import pdb
    pdb.set_trace()
    K, QS = kinship_process(K, nsamples, verbose)

    if lik == 'normal':
        model = _perform_lmm(
            _binomial_y(y.values, lik), M, QS, G, mixed, verbose)
    else:
        model = _perform_glmm(
            _binomial_y(y.values, lik), lik, M, K, QS, G, mixed, verbose)

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


def _perform_glmm(y, lik, M, K, QS, G, mixed, verbose):
    from pandas import Series

    glmm = GLMMExpFam(y, lik, M.values, QS)
    if not mixed or lik == 'bernoulli':
        glmm.delta = 1
        glmm.fix('delta')
    glmm.fit(verbose=verbose)

    eta = glmm.site.eta
    tau = glmm.site.tau

    gnormal = GLMMNormal(eta, tau, M.values, QS)
    gnormal.fit(verbose=verbose)

    beta = gnormal.beta

    keys = list(M.keys())
    ncov_effsizes = Series(beta, keys)

    flmm = gnormal.get_fast_scanner()
    flmm.set_scale(1.0)
    null_lml = flmm.null_lml()

    alt_lmls, effsizes = flmm.fast_scan(G.values, verbose=verbose)

    alt_lmls = Series(alt_lmls, list(G.keys()))
    effsizes = Series(effsizes, list(G.keys()))

    return QTLModel(null_lml, alt_lmls, effsizes, ncov_effsizes)


def _binomial_y(y, lik):
    # ugly hack, remove this when possible
    if lik == 'binomial':
        return y[:, 0], y[:, 1]
    return y.ravel()


def _intersect_indices(G, y, K, M):
    indices = G.index.copy()
    indices = indices.intersection(y.index)
    if K is not None:
        indices = indices.intersection(K.index)
    indices = indices.intersection(M.index)
    # return Index(list(indices))
    return indices
