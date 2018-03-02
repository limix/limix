from __future__ import division

from numpy import all as npall
from numpy import diag, isfinite, ones, asarray
from pandas import DataFrame

from glimix_core.glmm import GLMMExpFam, GLMMNormal
from glimix_core.lmm import LMM
from numpy_sugar.linalg import economic_qs

from .model import QTLModel
from ..nice_arrays import (assure_named_columns, covariates_process,
                           kinship_process, phenotype_process)

from ..nice_arrays import normalise_phenotype_matrix
from ..nice_arrays import normalise_covariates_matrix
from ..nice_arrays import normalise_kinship_matrix
from ..nice_arrays import normalise_candidates_matrix
from ..nice_arrays import infer_samples_index
from ..nice_arrays import default_covariates_index
from ..nice_arrays import default_candidates_index
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
        >>> print(model.variant_pvalues)  # doctest: +NPY_FLEX_NUMS
        rs0    0.554444
        rs1    0.218996
        rs2    0.552201
        dtype: float64
        >>> print(model.variant_effsizes)  # doctest: +NPY_FLEX_NUMS
        rs0   -0.130866
        rs1   -0.315077
        rs2   -0.143869
        dtype: float64
        >>> print(model.variant_effsizes_se)  # doctest: +NPY_FLEX_NUMS
        rs0    0.221389
        rs1    0.256327
        rs2    0.242013
        dtype: float64
        >>> print(model)  # doctest: +NPY_FLEX_NUMS
        Variants
               effsizes  effsizes_se   pvalues
        count  3.000000     3.000000  3.000000
        mean  -0.196604     0.239910  0.441880
        std    0.102807     0.017563  0.193027
        min   -0.315077     0.221390  0.218996
        25%   -0.229473     0.231701  0.385598
        50%   -0.143869     0.242013  0.552201
        75%   -0.137368     0.249170  0.553322
        max   -0.130866     0.256327  0.554444
        <BLANKLINE>
        Covariate effect sizes for the null model
              age   offset
        -0.005568  0.39529
    """
    lik = lik.lower()

    if verbose:
        print_analysis(lik, "Quantitative trait locus analysis")

    nsamples = len(G)
    if isinstance(y, (tuple, list)):
        y = asarray(y, float).T

    if M is None:
        M = ones((nsamples, 1))

    arrs = [y, G, M]
    if K is not None:
        arrs.append(K)

    samples_index = infer_samples_index(arrs)
    if len(samples_index) == 0:
        raise ValueError("Could not infer an index for samples."
                         " Please, check if the passed arrays are in order.")

    if hasattr(y, 'index'):
        y = y.loc[y.index.intersection(samples_index)]
    else:
        y = DataFrame(data=y, index=samples_index.copy())

    if hasattr(G, 'index'):
        G = G.loc[G.index.intersection(samples_index)]
    else:
        G = DataFrame(
            data=G,
            index=samples_index.copy(),
            columns=default_candidates_index(G.shape[1]))

    if hasattr(M, 'index'):
        M = M.loc[M.index.intersection(samples_index)]
    else:
        M = DataFrame(
            data=M,
            index=samples_index.copy(),
            columns=default_covariates_index(M.shape[1]))

    if K is not None:
        if hasattr(K, 'index'):
            K = K.loc[K.index.intersection(samples_index), :]
            K = K.loc[:, K.columns.intersection(samples_index)]
        else:
            K = DataFrame(
                data=K,
                index=samples_index.copy(),
                columns=samples_index.copy())

    y = normalise_phenotype_matrix(y, lik)
    if K is not None:
        K = normalise_kinship_matrix(K)

    G = normalise_candidates_matrix(G)
    M = covariates_process(M, nsamples)
    M = normalise_covariates_matrix(M)

    y = phenotype_process(lik, y)

    if not npall(isfinite(G)):
        raise ValueError("Variant values must be finite.")

    mixed = K is not None

    indices = _intersect_indices(G, y, K, M)

    G = G.loc[indices, :]
    y = y.loc[indices, :]

    if K is not None:
        K = K.loc[indices, :]
        K = K.loc[:, indices]

    M = M.loc[indices, :]

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
    print("GLMM: delta: {}, scale: {}".format(glmm.delta, glmm.scale))

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
