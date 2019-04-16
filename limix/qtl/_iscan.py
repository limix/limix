import sys

from limix._display import session_line

from .._bits import unvec
from .._data import asarray as _asarray, conform_dataset, normalize_likelihood
from .._display import session_block
from ._assert import assert_finite
from ._result import IScanResultFactory


def iscan(G, y, lik="normal", K=None, M=None, idx=None, E0=None, E1=None, verbose=True):
    r"""
    Single-trait association with interaction test via generalized linear mixed models.

    The general formulae for normally distributed traits is

    .. math::

        𝐲 = 𝙼𝛂 + (𝙶⊙𝙴₀)𝛃₀ + (𝙶⊙𝙴₁)𝛃₁ + 𝐮 + 𝛆,\\
        \text{where}~~ 𝐮∼𝓝(𝟎, 𝓋₀𝙺) ~~\text{and}~~ 𝛆∼𝓝(𝟎, 𝓋₁𝙸).

    The operator ⊙ works as follows:

    .. math::

        𝙰⊙𝙱 = [𝙰₀𝙱₀ ~~...~~ 𝙰₀𝙱ₙ ~~ 𝙰₁𝙱₀ ~~...~~ 𝙰₁𝙱ₙ ~~...~~ 𝙰ₘ𝙱ₙ]

    The covariates is enconded in matrix 𝙼 while the candidate set is enconded in matrix
    𝙶. The parameters are the effect sizes 𝛂, 𝛃₀, and 𝛃₁, and the variances 𝓋₀ and 𝓋₁.

    It performs likelihood-ratio tests for the following cases, where the first
    hypothesis is the null one while the second hypothesis is the alternative one:

    - H₀ vs H₁: testing for vec(𝛃₀) ≠ 𝟎 while vec(𝛃₁) = 𝟎
    - H₀ vs H₂: testing for [vec(𝛃₀) vec(𝛃₁)] ≠ 𝟎
    - H₁ vs H₂: testing for vec(𝛃₁) ≠ 𝟎

    It also supports generalized linear mixed models (GLMM). In this case, the following
    likelihoods are implemented:
    - Bernoulli
    - Probit
    - Binomial
    - Poisson

    Formally, let p(𝜇) be one of the supported probability distributions where 𝜇 is
    its mean. The H₀ model is defined as follows:

    .. math::

        yᵢ ∼ p(𝜇ᵢ=g(zᵢ)) ~~\text{for}~~ 𝐳 ∼ 𝓝(𝙼𝛂 + (𝙶⊙𝙴₀)𝛃₀ + (𝙶⊙𝙴₁)𝛃₁, 𝓋₀𝙺 + 𝓋₁𝙸).

    g(⋅) is the corresponding canonical link function for the Bernoulli, Binomial, and
    Poisson likelihoods. The Probit likelihood, on the other hand, is a Bernoulli
    likelihood with probit link function.

    Parameters
    ----------
    G : n×m array_like
        Genetic candidates.
    Y : n×p array_like
        Rows are samples and columns are phenotypes.
    lik : tuple, "normal", "bernoulli", "probit", "binomial", "poisson"
        Sample likelihood describing the residual distribution.
        Either a tuple or a string specifying the likelihood is required. The Normal,
        Bernoulli, Probit, and Poisson likelihoods can be selected by providing a
        string. Binomial likelihood on the other hand requires a tuple because of the
        number of trials: ``("binomial", array_like)``. Defaults to ``"normal"``.
    K : n×n array_like
        Sample covariance, often the so-called kinship matrix.
    M : n×c array_like
        Covariates matrix.
    idx : list
        List of candidate indices that defines the set of candidates to be used in the
        tests.
    E0 : array_like
        Matrix representing the first environment.
    E1 : array_like
        Matrix representing the second environment.
    verbose : bool, optional
        ``True`` to display progress and summary; ``False`` otherwise.

    Returns
    -------
    result : :class:`limix.qtl._result.IScanResult`
        P-values, log of marginal likelihoods, effect sizes, and associated statistics.

    Notes
    -----
    It will raise a ``ValueError`` exception if non-finite values are passed. Please,
    refer to the :func:`limix.qc.mean_impute` function for missing value imputation.
    """
    from numpy_sugar.linalg import economic_qs
    from xarray import concat
    from numpy import asarray, empty, ones

    lik = normalize_likelihood(lik)
    lik_name = lik[0]

    with session_block("QTL analysis", disable=not verbose):

        with session_line("Normalising input... ", disable=not verbose):

            data = conform_dataset(y, M, G=G, K=K)

        Y = data["y"]
        M = data["M"]
        G = data["G"]
        K = data["K"]

        assert_finite(y, M, K)
        nsamples = y.shape[0]

        if E1 is None:
            E1 = ones((nsamples, 1))

        if E0 is None:
            E0 = empty((nsamples, 0))

        E0 = _asarray(E0, "env0", ["sample", "env"])
        E1 = _asarray(E1, "env1", ["sample", "env"])
        E01 = concat([E0, E1], dim="env")

        if K is not None:
            QS = economic_qs(K)
        else:
            QS = None

        if lik_name == "normal":
            scanner, v0, v1 = _lmm(Y.values.ravel(), M.values, QS, verbose)
        else:
            scanner, v0, v1 = _glmm(Y.values.ravel(), lik, M.values, QS, verbose)

        r = IScanResultFactory(
            lik_name,
            Y.trait,
            M.covariate,
            G.candidate,
            E0.env,
            E1.env,
            scanner.null_lml,
            scanner.null_beta,
            scanner.null_beta_se,
            v0,
            v1,
        )

        if idx is None:

            assert E1.shape[1] > 0
            idx = range(G.shape[1])

            if E0.shape[1] == 0:
                r1 = scanner.fast_scan(G, verbose)

            for i in idx:
                i = _2d_sel(i)
                g = asarray(G[:, i], float)

                if E0.shape[1] > 0:
                    r1 = scanner.scan(g, E0)
                    h1 = _normalise_scan_names(r1)
                else:
                    h1 = _normalise_scan_names({k: v[i] for k, v in r1.items()})
                    h1["covariate_effsizes"] = h1["covariate_effsizes"].ravel()
                    h1["covariate_effsizes_se"] = h1["covariate_effsizes_se"].ravel()

                r2 = scanner.scan(g, E01)
                h2 = _normalise_scan_names(r2)
                r.add_test(i, h1, h2)
        else:
            for i in idx:
                i = _2d_sel(i)
                g = asarray(G[:, i], float)

                r1 = scanner.scan(g, E0)
                r2 = scanner.scan(g, E01)

                h1 = _normalise_scan_names(r1)
                h2 = _normalise_scan_names(r2)
                r.add_test(i, h1, h2)

        r = r.create()
        if verbose:
            print(r)

        return r


def _normalise_scan_names(r):
    return {
        "lml": r["lml"],
        "covariate_effsizes": r["effsizes0"],
        "covariate_effsizes_se": r["effsizes0_se"],
        "candidate_effsizes": r["effsizes1"],
        "candidate_effsizes_se": r["effsizes1_se"],
        "scale": r["scale"],
    }


def _2d_sel(idx):
    from collections.abc import Iterable

    if not isinstance(idx, (slice, Iterable)):
        return [idx]

    return idx


class ScannerWrapper:
    def __init__(self, scanner):
        self._scanner = scanner

    @property
    def null_lml(self):
        return self._scanner.null_lml()

    @property
    def null_beta(self):
        return self._scanner.null_beta

    @property
    def null_beta_se(self):
        from numpy import sqrt

        se = sqrt(self._scanner.null_beta_covariance.diagonal())
        return se

    def fast_scan(self, G, verbose):
        r = self._scanner.fast_scan(G, verbose=verbose)
        r["effsizes1"] = unvec(r["effsizes1"], (-1, G.shape[1])).T
        r["effsizes1_se"] = unvec(r["effsizes1_se"], (-1, G.shape[1])).T
        return r

    def scan(self, G, E):
        from .._bits import cdot

        r = self._scanner.scan(cdot(G, E))
        r["effsizes1"] = unvec(r["effsizes1"], (-1, G.shape[1])).T
        r["effsizes1_se"] = unvec(r["effsizes1_se"], (-1, G.shape[1])).T
        return r


def _lmm(y, M, QS, verbose):
    from glimix_core.lmm import LMM

    lmm = LMM(y, M, QS, restricted=False)
    lmm.fit(verbose=verbose)
    sys.stdout.flush()

    if QS is None:
        v0 = None
    else:
        v0 = lmm.v0
    v1 = lmm.v1
    scanner = ScannerWrapper(lmm.get_fast_scanner())

    return scanner, v0, v1


def _glmm(y, lik, M, QS, verbose):
    from glimix_core.glmm import GLMMExpFam, GLMMNormal

    glmm = GLMMExpFam(y.ravel(), lik, M, QS)

    glmm.fit(verbose=verbose)
    v0 = glmm.v0
    v1 = glmm.v1
    sys.stdout.flush()

    eta = glmm.site.eta
    tau = glmm.site.tau

    gnormal = GLMMNormal(eta, tau, M, QS)
    gnormal.fit(verbose=verbose)

    scanner = ScannerWrapper(gnormal.get_fast_scanner())

    return scanner, v0, v1
