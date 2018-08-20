from __future__ import division

import sys

from numpy import ones

from limix.display import timer_text

from .._dataset import normalise_dataset
from ..display import session_text
from .._likelihood import assert_likelihood_name, normalise_extreme_values
from ._model import QTLModel


def scan(G, y, lik, K=None, M=None, verbose=True):
    r"""Single-variant association testing via generalised linear mixed models.

    It supports Normal (linear mixed model), Bernoulli, Probit, Binomial, and Poisson
    residual errors, defined by ``lik``.
    The columns of ``G`` define the candidates to be tested for association
    with the phenotype ``y``.
    The covariance matrix is set by ``K``.
    If not provided, or set to ``None``, the generalised linear model
    without random effects is assumed.
    The covariates can be set via the parameter ``M``.
    We recommend to always provide a column of ones when covariates are actually
    provided.

    Parameters
    ----------
    G : array_like
        :math:`N` individuals by :math:`S` candidate markers.
    y : tuple, array_like
        Either a tuple of two arrays of :math:`N` individuals each (Binomial
        phenotypes) or an array of :math:`N` individuals (Normal, Poisson, Bernoulli, or
        Probit phenotypes).
    lik : "normal", "bernoulli", "probit", binomial", "poisson"
        Sample likelihood describing the residual distribution.
    K : array_like, optional
        :math:`N`-by-:math:`N` covariance matrix (e.g., kinship coefficients).
        Set to ``None`` for a generalised linear model without random effects.
        Defaults to ``None``.
    M : array_like, optional
        `N` individuals by `S` covariates.
        It will create a :math:`N`-by-:math:`1` matrix ``M`` of ones representing the
        offset covariate if ``None`` is passed. If an array is passed, it will used as
        is. Defaults to ``None``.
    verbose : bool, optional
        ``True`` to display progress and summary; ``False`` otherwise.

    Returns
    -------
    :class:`limix.qtl.QTLModel`
        QTL representation.

    Examples
    --------
    .. doctest::

        >>> from numpy import dot, exp, sqrt, ones
        >>> from numpy.random import RandomState
        >>> from pandas import DataFrame
        >>> import pandas as pd
        >>> from limix.qtl import scan
        >>>
        >>> random = RandomState(1)
        >>> pd.options.display.float_format = "{:9.6f}".format
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
        >>> model.variant_pvalues  # doctest: +FLOAT_CMP
        rs0    0.554443
        rs1    0.218996
        rs2    0.552200
        dtype: float64
        >>> model.variant_effsizes  # doctest: +FLOAT_CMP
        rs0   -0.130866
        rs1   -0.315077
        rs2   -0.143869
        dtype: float64
        >>> model.variant_effsizes_se  # doctest: +FLOAT_CMP
        rs0    0.221389
        rs1    0.256326
        rs2    0.242013
        dtype: float64
        >>> model  # doctest: +FLOAT_CMP
        Variants
        --------
               effsizes  effsizes_se   pvalues
        count  3.000000     3.000000  3.000000
        mean  -0.196604     0.239910  0.441880
        std    0.102807     0.017563  0.193027
        min   -0.315077     0.221389  0.218996
        25%   -0.229473     0.231701  0.385598
        50%   -0.143869     0.242013  0.552200
        75%   -0.137367     0.249170  0.553322
        max   -0.130866     0.256326  0.554443
        <BLANKLINE>
        Covariate effect sizes for H0
        -----------------------------
              age    offset
        -0.005568  0.395287

    Notes
    -----
    It will raise a ``ValueError`` exception if non-finite values are passed. Please,
    refer to the :func:`limix.qc.mean_impute` function for missing value imputation.
    """
    from numpy_sugar import is_all_finite
    from numpy_sugar.linalg import economic_qs

    if not isinstance(lik, (tuple, list)):
        lik = (lik,)

    lik_name = lik[0].lower()
    assert_likelihood_name(lik_name)

    if M is None:
        M = ones((len(y), 1))

    with session_text("qtl analysis", disable=not verbose):

        with timer_text("Normalising input...", disable=not verbose):
            data = normalise_dataset(y, M, G=G, K=K)

        y = data["y"]
        M = data["M"]
        G = data["G"]
        K = data["K"]

        if not is_all_finite(y):
            raise ValueError("Outcome must have finite values only.")

        if not is_all_finite(M):
            raise ValueError("Covariates must have finite values only.")

        if K is not None:
            if not is_all_finite(K):
                raise ValueError("Covariate matrix must have finite values only.")
            QS = economic_qs(K)
        else:
            QS = None

        y = normalise_extreme_values(data["y"], lik)

        if lik_name == "normal":
            model = _perform_lmm(y.values, M, QS, G, verbose)
        else:
            model = _perform_glmm(y.values, lik, M, K, QS, G, verbose)

        return model


def _perform_lmm(y, M, QS, G, verbose):
    from glimix_core.lmm import LMM
    from pandas import Series

    lmm = LMM(y, M.values, QS)

    lmm.fit(verbose=verbose)
    sys.stdout.flush()

    null_lml = lmm.lml()

    beta = lmm.beta

    covariates = list(M.coords["covariate"].values)
    ncov_effsizes = Series(beta, covariates)

    flmm = lmm.get_fast_scanner()
    alt_lmls, effsizes = flmm.fast_scan(G.values, verbose=verbose)

    candidates = list(G.coords["candidate"].values)
    alt_lmls = Series(alt_lmls, candidates)
    effsizes = Series(effsizes, candidates)

    return QTLModel(null_lml, alt_lmls, effsizes, ncov_effsizes)


def _perform_glmm(y, lik, M, K, QS, G, verbose):
    from glimix_core.glmm import GLMMExpFam, GLMMNormal
    from pandas import Series, DataFrame

    glmm = GLMMExpFam(y.ravel(), lik, M.values, QS)
    glmm.fit(verbose=verbose)
    sys.stdout.flush()

    eta = glmm.site.eta
    tau = glmm.site.tau

    gnormal = GLMMNormal(eta, tau, M.values, QS)
    gnormal.fit(verbose=verbose)

    beta = gnormal.beta

    covariates = list(M.coords["covariate"].values)
    ncov_effsizes = Series(beta, covariates)

    flmm = gnormal.get_fast_scanner()
    flmm.set_scale(1.0)
    null_lml = flmm.null_lml()

    if isinstance(G, DataFrame):
        alt_lmls, effsizes = flmm.fast_scan(G.values, verbose=verbose)
    else:
        alt_lmls, effsizes = flmm.fast_scan(G, verbose=verbose)

    candidates = list(G.coords["candidate"].values)
    alt_lmls = Series(alt_lmls, candidates)
    effsizes = Series(effsizes, candidates)

    return QTLModel(null_lml, alt_lmls, effsizes, ncov_effsizes)
