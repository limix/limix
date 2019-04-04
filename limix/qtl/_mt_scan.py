import sys

from limix._display import session_line

from .._data import asarray as _asarray, conform_dataset
from .._display import session_block
from ._result import ScanResultFactory


def mt_scan(G, Y, idx=None, K=None, M=None, A=None, A0=None, A1=None, verbose=True):
    """
    Multi-trait association and interaction testing via linear mixed models.

    Let n, c, and p be the number of samples, covariates, and traits, respectively.
    The outcome variable Y is a n×p matrix distributed according to ::

        vec(Y) ~ N((A ⊗ M) vec(B), K₀ = C₀ ⊗ K + C₁ ⊗ I) under H₀.

    A and M are design matrices of dimensions p×p and n×c provided by the user,
    where X is the usual matrix of covariates commonly used in single-trait models.
    B is a c×p matrix of fixed-effect sizes per trait.
    C₀ and C₁ are both symmetric matrices of dimensions p×p, for which C₁ is
    guaranteed by our implementation to be of full rank.
    The parameters of the H₀ model are the matrices B, C₀, and C₁.

    The additional models H₁ and H₂ are define as ::

        vec(Y) ~ N((A ⊗ M) vec(B) + (A₀ ⊗ Gᵢ) vec(A₁), s⋅K₀)

    and ::

        vec(Y) ~ N((A ⊗ M) vec(B) + (A₀ ⊗ Gᵢ) vec(A₁) + (A₁ ⊗ Gᵢ) vec(A₂), s⋅K₀)

    It performs likelihood-ratio tests for the following cases, where the first
    hypothesis is the null one while the second hypothesis is the alternative one:
    - H₀ vs H₁: testing for vec(A₁) ≠ 0 while vec(A₂) = 0
    - H₀ vs H₂: testing for [vec(A₁) vec(A₂)] ≠ 0
    - H₁ vs H₂: testing for vec(A₂) ≠ 0

    Parameters
    ----------
    G : n×m array_like
        Genetic candidates.
    Y : n×p array_like
        p phenotype values for n samples.
    idx : list
        List of candidate indices that defines the set of candidates to be used in the
        tests.
    K : n×n array_like
        Sample covariance, often the so-called kinship matrix.
    M : n×c array_like
        Covariates matrix.
    A : p×p array_like
        Symmetric trait-by-trait design matrix.
    A0 : p×p₀ array_like
        Matrix A₀, possibility a non-symmetric one.
    A1 : p×p₁ array_like
        Matrix A₁, possibility a non-symmetric one.
    verbose : bool, optional
        ``True`` to display progress and summary; ``False`` otherwise.

    Returns
    -------
    :class:`limix.qtl.ScanResult`
        P-values, log of marginal likelihoods, effect sizes, and associated statistics.

    Examples
    --------
    .. doctest::

        >>> from limix.qtl import mt_scan
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
        >>> r = mt_scan(G, Y, idx=idx, K=K, M=M, A=A, A0=A0, A1=A1, verbose=False)
    """
    from glimix_core.lmm import Kron2Sum
    from numpy_sugar.linalg import economic_qs, ddot
    from xarray import concat
    from ._assert import assert_finite
    from numpy import eye, sqrt, asarray, empty

    with session_block("qtl analysis", disable=not verbose):

        with session_line("Normalising input... ", disable=not verbose):

            data = conform_dataset(Y, M, G=G, K=K)

        Y = data["y"]
        M = data["M"]
        G = data["G"]
        K = data["K"]

        assert_finite(Y, M, K)

        if A is None:
            A = eye(Y.shape[1])

        if A1 is None:
            A1 = eye(Y.shape[1])

        if A0 is None:
            A0 = empty((Y.shape[1], 0))

        A0 = _asarray(A0, "env0", ["sample", "env"])
        A1 = _asarray(A1, "env1", ["sample", "env"])
        A01 = concat([A0, A1], dim="env")

        if idx is None:
            idx = range(G.shape[1])

        if K is not None:
            QS = economic_qs(K)
        else:
            QS = None

        KG = ddot(QS[0][0], sqrt(QS[1]))
        lmm = Kron2Sum(Y.values, A, M.values, KG, restricted=False)
        lmm.fit(verbose=verbose)
        sys.stdout.flush()

        scanner = lmm.get_fast_scanner()
        r = ScanResultFactory(
            "normal",
            Y.trait,
            M.covariate,
            G.candidate,
            A0.env,
            A1.env,
            scanner.null_lml(),
            scanner.null_beta,
            sqrt(scanner.null_beta_covariance.diagonal()),
            lmm.C0,
            lmm.C1,
        )

        for i in idx:

            i = _2d_sel(i)
            g = asarray(G[:, i], float)

            r1 = scanner.scan(A0, g)
            r2 = scanner.scan(A01, g)

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
