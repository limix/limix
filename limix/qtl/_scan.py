import sys

from limix._display import session_line

from .._data import asarray as _asarray, conform_dataset, normalize_likelihood
from .._display import session_block
from ._assert import assert_finite
from ._result import MTScanResultFactory, STScanResultFactory


def scan(
    G, Y, lik="normal", K=None, M=None, idx=None, A=None, A0=None, A1=None, verbose=True
):
    """
    Multi-trait association and interaction testing via linear mixed models.

    Let n, c, and p be the number of samples, covariates, and traits, respectively.
    The outcome variable Y is a n×p matrix distributed according to ::

        vec(Y) ~ N((A ⊗ M) vec(𝚨), K₀ = C₀ ⊗ K + C₁ ⊗ I) under H₀.

    A and M are design matrices of dimensions p×p and n×c provided by the user,
    where X is the usual matrix of covariates commonly used in single-trait models.
    𝚨 is a c×p matrix of fixed-effect sizes per trait.
    C₀ and C₁ are both symmetric matrices of dimensions p×p, for which C₁ is
    guaranteed by our implementation to be of full rank.
    The parameters of the H₀ model are the matrices 𝚨, C₀, and C₁.

    The additional models H₁ and H₂ are define as ::

        vec(Y) ~ N((A ⊗ M) vec(𝚨) + (A₀ ⊗ Gᵢ) vec(𝚩₀), s⋅K₀)

    and ::

        vec(Y) ~ N((A ⊗ M) vec(𝚨) + (A₀ ⊗ Gᵢ) vec(𝚩₀) + (A₁ ⊗ Gᵢ) vec(𝚩₁), s⋅K₀)

    It performs likelihood-ratio tests for the following cases, where the first
    hypothesis is the null one while the second hypothesis is the alternative one:

    - H₀ vs H₁: testing for vec(𝚩₀) ≠ 𝟎 while vec(𝚩₁) = 𝟎
    - H₀ vs H₂: testing for [vec(𝚩₀) vec(𝚩₁)] ≠ 𝟎
    - H₁ vs H₂: testing for vec(𝚩₁) ≠ 𝟎

    It supports generalized linear mixed models (GLMM) when a single trait is used.
    In this case, the following likelihoods are implemented:

    - Bernoulli
    - Probit
    - Binomial
    - Poisson

    Formally, let p(𝜇) be one of the supported probability distributions where 𝜇 is
    its mean. The H₀ model is defined as follows::

        yᵢ ∼ p(𝜇ᵢ=g(zᵢ)) for 𝐳 ∼ 𝓝(..., ...).

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
    A : p×p array_like
        Symmetric trait-by-trait design matrix.
    A0 : p×p₀ array_like, optional
        Matrix A₀, possibility a non-symmetric one. If ``None``, it defines an empty
        matrix, p₀=0. Defaults to ``None``.
    A1 : p×p₁ array_like, optional
        Matrix A₁, possibility a non-symmetric one. If ``None``, it defines an identity
        matrix, p₀=p. Defaults to ``None``.
    verbose : bool, optional
        ``True`` to display progress and summary; ``False`` otherwise.

    Returns
    -------
    result : :class:`limix.qtl._result.STScanResult`, :class:`limix.qtl._result.MTScanResult`
        P-values, log of marginal likelihoods, effect sizes, and associated statistics.

    Examples
    --------
    .. doctest::

        >>> from limix.qtl import scan
        >>> from numpy import reshape, kron, eye
        >>> from numpy import concatenate
        >>> from numpy.random import RandomState
        >>> import scipy.stats as st
        >>> from limix.qc import normalise_covariance
        >>>
        >>> def vec(x):
        ...     return reshape(x, (-1,) + x.shape[2:], order="F")
        >>>
        >>> def unvec(x, shape):
        ...     return reshape(x, shape, order="F")
        >>>
        >>> random = RandomState(0)
        >>> n = 30
        >>> ntraits = 2
        >>> ncovariates = 3
        >>>
        >>> A = random.randn(ntraits, ntraits)
        >>> A = A @ A.T
        >>> M = random.randn(n, ncovariates)
        >>>
        >>> C0 = random.randn(ntraits, ntraits)
        >>> C0 = C0 @ C0.T
        >>>
        >>> C1 = random.randn(ntraits, ntraits)
        >>> C1 = C1 @ C1.T
        >>>
        >>> G = random.randn(n, 4)
        >>>
        >>> A0 = random.randn(ntraits, 1)
        >>> A1 = random.randn(ntraits, 2)
        >>> A01 = concatenate((A0, A1), axis=1)
        >>>
        >>> K = random.randn(n, n + 1)
        >>> K = normalise_covariance(K @ K.T)
        >>>
        >>> beta = vec(random.randn(ntraits, ncovariates))
        >>> alpha = vec(random.randn(A01.shape[1], G.shape[1]))
        >>>
        >>> mvn = st.multivariate_normal
        >>> m = kron(A, M) @ beta + kron(A01, G) @ alpha
        >>> Y = unvec(mvn(m, kron(C0, K) + kron(C1, eye(n))).rvs(), (n, -1))
        >>>
        >>> idx = [[0, 1], 2, [3]]
        >>> r = scan(G, Y, idx=idx, K=K, M=M, A=A, A0=A0, A1=A1, verbose=False)

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
        >>> result = scan(candidates, y, 'poisson', K, M=M, verbose=False)
        >>>
        >>> result.stats  # doctest: +FLOAT_CMP
                   lml0       lml2  dof20    scale2      pv20
        test
        0    -48.720890 -48.536860      1  0.943532  0.544063
        1    -48.720890 -47.908341      1  0.904814  0.202382
        2    -48.720890 -48.534754      1  0.943400  0.541768
        >>> print(result)  # doctest: +FLOAT_CMP
        Hypothesis 0
        ------------
        <BLANKLINE>
        𝐳 ~ 𝓝(𝙼𝜶, 0.000⋅𝙺 + 0.788⋅𝙸) for yᵢ ~ Poisson(λᵢ=g(zᵢ)) and g(x)=eˣ
        <BLANKLINE>
        M     = ['offset' 'age']
        𝜶     = [ 0.39528889 -0.00556797]
        se(𝜶) = [0.50173695 0.01505240]
        lml   = -48.720890273519444
        <BLANKLINE>
        Hypothesis 2
        ------------
        <BLANKLINE>
        𝐳 ~ 𝓝(𝙼𝜶 + G𝛃, s(0.000⋅𝙺 + 0.788⋅𝙸)) for yᵢ ~ Poisson(λᵢ=g(zᵢ)) and g(x)=eˣ
        <BLANKLINE>
                  lml       cov. effsizes   cand. effsizes
        --------------------------------------------------
        mean   -4.833e+01       2.393e-01       -1.966e-01
        std     3.623e-01       2.713e-01        1.028e-01
        min    -4.854e+01      -8.490e-03       -3.151e-01
        25%    -4.854e+01      -7.684e-03       -2.295e-01
        50%    -4.853e+01       2.243e-01       -1.439e-01
        75%    -4.822e+01       4.725e-01       -1.374e-01
        max    -4.791e+01       5.255e-01       -1.309e-01
        <BLANKLINE>
        Likelihood-ratio test p-values
        ------------------------------
        <BLANKLINE>
               𝓗₀ vs 𝓗₂
        ----------------
        mean   4.294e-01
        std    1.966e-01
        min    2.024e-01
        25%    3.721e-01
        50%    5.418e-01
        75%    5.429e-01
        max    5.441e-01
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
        >>> result = scan(X, successes, ("binomial", ntrials), K, verbose=False)
        >>> print(result)  # doctest: +FLOAT_CMP
        Hypothesis 0
        ------------
        <BLANKLINE>
        𝐳 ~ 𝓝(𝙼𝜶, 0.152⋅𝙺 + 1.738⋅𝙸) for yᵢ ~ Binom(μᵢ=g(zᵢ), nᵢ) and g(x)=1/(1+e⁻ˣ)
        <BLANKLINE>
        M     = ['offset']
        𝜶     = [0.40956942]
        se(𝜶) = [0.55141166]
        lml   = -142.80784719977515
        <BLANKLINE>
        Hypothesis 2
        ------------
        <BLANKLINE>
        𝐳 ~ 𝓝(𝙼𝜶 + G𝛃, s(0.152⋅𝙺 + 1.738⋅𝙸)) for yᵢ ~ Binom(μᵢ=g(zᵢ), nᵢ) and g(x)=1/(1+e⁻ˣ)
        <BLANKLINE>
                  lml       cov. effsizes   cand. effsizes
        --------------------------------------------------
        mean   -1.425e+02       3.701e-01        2.271e-01
        std     4.110e-01       2.296e-02        5.680e-01
        min    -1.427e+02       3.539e-01       -1.745e-01
        25%    -1.426e+02       3.620e-01        2.631e-02
        50%    -1.425e+02       3.701e-01        2.271e-01
        75%    -1.423e+02       3.782e-01        4.279e-01
        max    -1.422e+02       3.864e-01        6.287e-01
        <BLANKLINE>
        Likelihood-ratio test p-values
        ------------------------------
        <BLANKLINE>
               𝓗₀ vs 𝓗₂
        ----------------
        mean   4.959e-01
        std    3.362e-01
        min    2.582e-01
        25%    3.771e-01
        50%    4.959e-01
        75%    6.148e-01
        max    7.336e-01

    Notes
    -----
    It will raise a ``ValueError`` exception if non-finite values are passed. Please,
    refer to the :func:`limix.qc.mean_impute` function for missing value imputation.
    """
    from numpy_sugar.linalg import economic_qs

    lik = normalize_likelihood(lik)

    if A is None:
        if A0 is not None or A1 is not None:
            raise ValueError("You cannot define `A0` or `A1` without defining `A`.")

    with session_block("QTL analysis", disable=not verbose):

        with session_line("Normalising input... ", disable=not verbose):
            data = conform_dataset(Y, M, G=G, K=K)

        Y = data["y"]
        M = data["M"]
        G = data["G"]
        K = data["K"]

        assert_finite(Y, M, K)

        if K is not None:
            QS = economic_qs(K)
        else:
            QS = None

        if verbose:
            print()
            _print_input_info(idx, lik, Y, M, G, K)
            print()

        if A is None:
            r = _single_trait_scan(idx, lik, Y, M, G, QS, verbose)
        else:
            r = _multi_trait_scan(idx, lik, Y, M, G, QS, A, A0, A1, verbose)

        r = r.create()
        if verbose:
            print()
            print(r)

        return r


def _print_input_info(idx, lik, Y, M, G, K):
    from limix._display import draw_list
    from limix._display import AlignedText, draw_title

    aligned = AlignedText(": ")
    likname = lik[0]
    aligned.add_item("Likelihood", likname)
    ntraits = Y.shape[1]
    traits = draw_list(Y.trait.values.tolist(), 5)
    aligned.add_item(f"Traits ({ntraits})", traits)

    ncovariates = M.shape[1]
    covariates = draw_list(M.covariate.values.tolist(), 5)
    aligned.add_item(f"Covariates ({ncovariates})", covariates)

    nvariants = G.shape[1]
    variants = draw_list(G.candidate.values.tolist(), 5)
    aligned.add_item(f"Variants {nvariants}", variants)
    if idx is None:
        ncandidates = nvariants
    else:
        ncandidates = len(idx)
    aligned.add_item("N. of candidates", ncandidates)

    if K is None:
        kinship_presence = "absent"
    else:
        kinship_presence = "present"
    aligned.add_item("Kinship", kinship_presence)

    print(draw_title("Input"))
    print(aligned.draw())


def _single_trait_scan(idx, lik, Y, M, G, QS, verbose):
    from numpy import asarray
    from tqdm import tqdm

    if lik[0] == "normal":
        scanner, v0, v1 = _st_lmm(Y.values.ravel(), M.values, QS, verbose)
    else:
        scanner, v0, v1 = _st_glmm(Y.values.ravel(), lik, M.values, QS, verbose)
    pass

    r = STScanResultFactory(
        lik[0],
        Y.trait.item(),
        M.covariate,
        G.candidate,
        scanner.null_lml(),
        scanner.null_beta,
        scanner.null_beta_se,
        v0,
        v1,
    )

    if idx is None:
        r1 = scanner.fast_scan(G, verbose)
        for i in tqdm(range(G.shape[1]), "Results", disable=not verbose):
            h2 = _normalise_scan_names({k: v[i] for k, v in r1.items()})
            r.add_test(i, h2)
    else:
        for i in tqdm(idx, "Results", disable=not verbose):
            i = _2d_sel(i)
            h2 = _normalise_scan_names(scanner.scan(asarray(G[:, i], float)))
            r.add_test(i, h2)
    return r


def _multi_trait_scan(idx, lik, Y, M, G, QS, A, A0, A1, verbose):
    from xarray import concat, DataArray
    from numpy import eye, asarray, empty
    from tqdm import tqdm

    ntraits = Y.shape[1]

    if A1 is None:
        A1 = eye(ntraits)
        A1 = DataArray(A1, dims=["sample", "env"], coords={"env": Y.trait.values})

    if A0 is None:
        A0 = empty((ntraits, 0))
        A0 = DataArray(A0, dims=["sample", "env"], coords={"env": asarray([], str)})

    A0 = _asarray(A0, "env0", ["sample", "env"])
    if "env" not in A0.coords:
        A0.coords["env"] = [f"env0_{i}" for i in range(A0.shape[1])]

    A1 = _asarray(A1, "env1", ["sample", "env"])
    if "env" not in A1.coords:
        A1.coords["env"] = [f"env1_{i}" for i in range(A1.shape[1])]

    A01 = concat([A0, A1], dim="env")

    if lik[0] == "normal":
        scanner, C0, C1 = _mt_lmm(Y, A, M, QS, verbose)
    else:
        msg = "Non-normal likelihood inference has not been implemented for"
        msg += " multiple traits yet."
        raise ValueError(msg)

    r = MTScanResultFactory(
        lik[0],
        Y.trait,
        M.covariate,
        G.candidate,
        A0.env,
        A1.env,
        scanner.null_lml(),
        scanner.null_beta,
        scanner.null_beta_se,
        C0,
        C1,
    )

    if idx is None:
        idx = range(G.shape[1])

    for i in tqdm(idx, "Results", disable=not verbose):

        i = _2d_sel(i)
        g = asarray(G[:, i], float)

        if A0.shape[1] == 0:
            h1 = None
        else:
            h1 = _normalise_scan_names(scanner.scan(A0, g))

        h2 = _normalise_scan_names(scanner.scan(A01, g))
        r.add_test(i, h1, h2)

    return r


def _st_lmm(Y, M, QS, verbose):
    from numpy import nan
    from glimix_core.lmm import LMM

    lmm = LMM(Y, M, QS, restricted=False)
    lmm.fit(verbose=verbose)
    sys.stdout.flush()

    if QS is None:
        v0 = nan
    else:
        v0 = lmm.v0

    v1 = lmm.v1

    return lmm.get_fast_scanner(), v0, v1


def _st_glmm(y, lik, M, QS, verbose):
    from numpy import nan
    from glimix_core.glmm import GLMMExpFam, GLMMNormal

    glmm = GLMMExpFam(y, lik, M, QS)

    glmm.fit(verbose=verbose)

    if QS is None:
        v0 = nan
    else:
        v0 = glmm.v0

    v1 = glmm.v1
    sys.stdout.flush()

    eta = glmm.site.eta
    tau = glmm.site.tau

    gnormal = GLMMNormal(eta, tau, M, QS)
    gnormal.fit(verbose=verbose)

    return gnormal.get_fast_scanner(), v0, v1


def _mt_lmm(Y, A, M, QS, verbose):
    from glimix_core.lmm import Kron2Sum
    from numpy_sugar.linalg import ddot
    from numpy import sqrt, zeros

    if QS is None:
        KG = zeros((Y.shape[0], 1))
    else:
        KG = ddot(QS[0][0], sqrt(QS[1]))

    lmm = Kron2Sum(Y.values, A, M.values, KG, restricted=False)
    lmm.fit(verbose=verbose)
    sys.stdout.flush()

    C0 = lmm.C0
    C1 = lmm.C1

    return lmm.get_fast_scanner(), C0, C1


def _2d_sel(idx):
    from collections.abc import Iterable

    if not isinstance(idx, (slice, Iterable)):
        return [idx]

    return idx


def _normalise_scan_names(r):
    from ._result._tuples import VariantResult

    return VariantResult(
        lml=r["lml"],
        covariate_effsizes=r["effsizes0"],
        candidate_effsizes=r["effsizes1"],
        covariate_effsizes_se=r["effsizes0_se"],
        candidate_effsizes_se=r["effsizes1_se"],
        scale=r["scale"],
    )
