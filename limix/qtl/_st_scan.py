from __future__ import division

import sys

from limix._display import session_line

from .._data import conform_dataset
from .._display import session_block
from .._data import check_likelihood_name
from ..qc._lik import normalise_extreme_values
from ._model import QTLModel


def st_scan(G, y, lik, K=None, M=None, verbose=True):
    r""" Single-variant association testing via generalised linear mixed models.

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
    y : array_like
        An outcome array of :math:`N` individuals.
    lik : tuple, "normal", "bernoulli", "probit", binomial", "poisson"
        Sample likelihood describing the residual distribution.
        Either a tuple or a string specifiying the likelihood is required. The Normal,
        Bernoulli, Probit, and Poisson likelihoods can be selected by providing a
        string. Binomial likelihood on the other hand requires a tuple because of the
        number of trials: ``("binomial", array_like)``.
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
        >>> from limix.qtl import st_scan
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
        >>> model = st_scan(candidates, y, 'poisson', K, M=M, verbose=False)
        >>>
        >>> model.variant_pvalues.to_dataframe()  # doctest: +FLOAT_CMP
                         pv
        candidate
        rs0        0.554444
        rs1        0.218996
        rs2        0.552200
        >>> model.variant_effsizes.to_dataframe()  # doctest: +FLOAT_CMP
                   effsizes
        candidate
        rs0       -0.130867
        rs1       -0.315078
        rs2       -0.143869
        >>> model.variant_effsizes_se.to_dataframe()  # doctest: +FLOAT_CMP
                   effsizes std
        candidate
        rs0            0.221390
        rs1            0.256327
        rs2            0.242013
        >>> model  # doctest: +FLOAT_CMP
        Variants
        --------
               effsizes  effsizes_se   pvalues
        count         3            3         3
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

    >>> from numpy import zeros
    >>>
    >>> nsamples = 50
    >>>
    >>> X = random.randn(nsamples, 2)
    >>> G = random.randn(nsamples, 100)
    >>> K = dot(G, G.T)
    >>> ntrials = random.randint(1, 100, nsamples)
    >>> z = dot(G, random.randn(100)) / sqrt(100)
    >>>
    >>> successes = zeros(len(ntrials), int)
    >>> for i, nt in enumerate(ntrials):
    ...     for _ in range(nt):
    ...         successes[i] += int(z[i] + 0.5 * random.randn() > 0)
    >>>
    >>> result = st_scan(X, successes, ("binomial", ntrials), K, verbose=False)
    >>> print(result)  # doctest: +FLOAT_CMP
    Variants
    --------
           effsizes  effsizes_se   pvalues
    count         2            2         2
    mean   0.227116     0.509575  0.478677
    std    0.567975     0.031268  0.341791
    min   -0.174503     0.487466  0.236994
    25%    0.026307     0.498520  0.357835
    50%    0.227116     0.509575  0.478677
    75%    0.427925     0.520630  0.599518
    max    0.628735     0.531685  0.720359
    <BLANKLINE>
    Covariate effect sizes for H0
    -----------------------------
       offset
     0.409570


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
    lik = (lik_name,) + lik[1:]
    check_likelihood_name(lik_name)

    with session_block("qtl analysis", disable=not verbose):

        with session_line("Normalising input... ", disable=not verbose):
            data = conform_dataset(y, M, G=G, K=K)

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

        if verbose:
            print(model)

        return model


def _perform_lmm(y, M, QS, G, verbose):
    from glimix_core.lmm import LMM
    from pandas import Series
    from xarray import DataArray

    lmm = LMM(y, M.values, QS)

    lmm.fit(verbose=verbose)
    sys.stdout.flush()

    null_lml = lmm.lml()

    beta = lmm.beta

    covariates = list(M.coords["covariate"].values)
    ncov_effsizes = Series(beta, covariates)

    flmm = lmm.get_fast_scanner()
    if hasattr(G, "data"):
        values = G.data
    else:
        values = G.values
    alt_lmls, effsizes = flmm.fast_scan(values, verbose=verbose)

    coords = {
        k: ("candidate", G.coords[k].values)
        for k in G.coords.keys()
        if G.coords[k].dims[0] == "candidate"
    }

    alt_lmls = DataArray(alt_lmls, dims=["candidate"], coords=coords)
    effsizes = DataArray(effsizes, dims=["candidate"], coords=coords)

    return QTLModel(null_lml, alt_lmls, effsizes, ncov_effsizes)


def _perform_glmm(y, lik, M, K, QS, G, verbose):
    from glimix_core.glmm import GLMMExpFam, GLMMNormal
    from pandas import Series
    from xarray import DataArray

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

    if hasattr(G, "data"):
        values = G.data
    else:
        values = G.values
    alt_lmls, effsizes = flmm.fast_scan(values, verbose=verbose)

    coords = {
        k: ("candidate", G.coords[k].values)
        for k in G.coords.keys()
        if G.coords[k].dims[0] == "candidate"
    }

    alt_lmls = DataArray(alt_lmls, dims=["candidate"], coords=coords)
    effsizes = DataArray(effsizes, dims=["candidate"], coords=coords)

    return QTLModel(null_lml, alt_lmls, effsizes, ncov_effsizes)
