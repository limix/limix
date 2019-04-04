import sys

from .._display import session_line
from .._data import conform_dataset
from .._display import session_block
from ._assert import assert_finite
from .._data import assert_likelihood
from ..qc._lik import normalise_extreme_values
from ._result import ScanResultFactory


def st_scan(
    G, y, lik, idx=None, K=None, M=None, A=None, A0=None, A1=None, verbose=True
):
    """
    Single-variant association testing via generalised linear mixed models.

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
    """
    from numpy_sugar.linalg import economic_qs

    if not isinstance(lik, (tuple, list)):
        lik = (lik,)

    lik_name = lik[0].lower()
    lik = (lik_name,) + lik[1:]
    assert_likelihood(lik_name)

    with session_block("qtl analysis", disable=not verbose):

        with session_line("Normalising input... ", disable=not verbose):
            data = conform_dataset(y, M, G=G, K=K)

        y = data["y"]
        M = data["M"]
        G = data["G"]
        K = data["K"]

        assert_finite(y, M, K)

        if K is not None:
            QS = economic_qs(K)
        else:
            QS = None

        y = normalise_extreme_values(data["y"], lik)

        if lik_name == "normal":
            r = _perform_lmm(y, idx, M, QS, G, verbose)
        else:
            r = _perform_glmm(y, lik, idx, M, K, QS, G, verbose)

        if verbose:
            print(r)

        return r


def _perform_lmm(y, idx, M, QS, G, verbose):
    from glimix_core.lmm import LMM
    from numpy import sqrt

    lmm = LMM(y.values, M.values, QS)

    lmm.fit(verbose=verbose)
    sys.stdout.flush()

    scanner = lmm.get_fast_scanner()

    r = ScanResultFactory(
        "normal",
        y.trait,
        M.covariate,
        G.candidate,
        [],
        ["0"],
        scanner.null_lml(),
        scanner.null_beta,
        sqrt(scanner.null_beta_covariance.diagonal()),
        lmm.v0,
        lmm.v1,
    )

    null_lml = lmm.lml()
    beta = lmm.beta

    # r.set_null(null_lml, beta, lmm.v1, lmm.v0)
    alt_lmls, effsizes = scanner.fast_scan(G.data, verbose=verbose)

    for i, data in enumerate(zip(alt_lmls, effsizes)):
        r.add_test(i, data[1], data[0])

    return r.create()


def _perform_glmm(y, lik, idx, M, K, QS, G, verbose):
    from glimix_core.glmm import GLMMExpFam, GLMMNormal

    glmm = GLMMExpFam(y.ravel(), lik, M.values, QS)

    glmm.fit(verbose=verbose)
    sys.stdout.flush()

    eta = glmm.site.eta
    tau = glmm.site.tau

    gnormal = GLMMNormal(eta, tau, M.values, QS)
    gnormal.fit(verbose=verbose)

    beta = gnormal.beta

    scanner = gnormal.get_fast_scanner()
    scanner.set_scale(1.0)
    null_lml = scanner.null_lml()

    r.set_null(null_lml, beta, gnormal.v1, gnormal.v0)

    alt_lmls, effsizes = scanner.fast_scan(G.values, verbose=verbose)

    for i, data in enumerate(zip(alt_lmls, effsizes)):
        r.add_test(i, data[1], data[0])

    return r.create()
