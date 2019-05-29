from limix._display import session_line

from .._data import conform_dataset, asarray as _asarray
from .._display import session_block
from ._assert import assert_finite


def sscan(G, y, E, M=None, idx=None, tests=None, verbose=True):
    r"""
    Structured linear mixed model that accounts for genotype-environment interactions.

    StructLMM [MC18]_ extends the conventional linear mixed model by including an
    additional per-individual effect term that accounts for genotype-environment
    interaction, which can be represented as an n×1 vector, 𝛃₁.
    The model can be cast as

    .. math::

        𝐲 = 𝙼𝛂 + 𝐠𝛽₀ + 𝐠⊙𝛃₁ + 𝐞 + 𝛆, ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
        \text{where}~~ 𝛃₁∼𝓝(𝟎, 𝓋₀Σ),~~ 𝐞∼𝓝(𝟎, 𝓋₁Σ),~~\text{and}~~ 𝛆∼𝓝(𝟎, 𝓋₂𝙸),

    where 𝛽₀ denotes the effect size of a conventional persistent genetic effect
    component and 𝛃₁ is a vector of per-individual effect sizes to account for
    heterogeneous genetic effects. The same environmental covariance Σ is used to
    account for genetic-environment interaction and for additive environmental effects.
    Vector 𝐠 is the variant of interest to test for association and/or interaction.
    The term 𝐠𝛽₀ represents the genetic effect, while term 𝐠⊙𝛃₁ represents the
    genetic-environment interaction effect.
    The parameters of the model are 𝛽₀, 𝛃₁, 𝓋₀, 𝓋₁, and 𝓋₂.

    It performs score tests for association and interaction, respectively::

    - H₀ vs H₁: testing for 𝛽₀ ≠ 0 while 𝛃₁ = 𝟎
    - H₀ vs H₂: testing for [𝛽₀ 𝛃₁] ≠ 𝟎

    Parameters
    ----------
    G : n×m array_like
        Genetic candidates.
    Y : n×p array_like
        Rows are samples and columns are phenotypes.
    E : n×𝚔 array_like
        Samples-by-environments design matrix.
    M : n×c array_like, optional
        Covariates matrix.
    idx : list, optional
        List of candidate indices that defines the set of candidates to be used in the
        tests.
        It defaults to ``range(n)`` such that every column of matrix G
    tests : list, optional
        List of tests to be performed.
        The possible values are ``"inter"`` and ``"assoc"``.
        It defaults to ``["inter"]``.
    verbose : bool, optional
        ``True`` to display progress and summary; ``False`` otherwise.

    References
    ----------
    .. [MC18] Moore, R., Casale, F. P., Bonder, M. J., Horta, D., Franke, L., Barroso,
        I., & Stegle, O. (2018). A linear mixed-model approach to study multivariate
        gene–environment interactions (p. 1). Nature Publishing Group.
    """
    from struct_lmm import StructLMM
    from numpy import asarray, zeros, hstack

    if tests is None:
        tests = set(["inter"])
    else:
        tests = set(tests)

    remain = tests - set(["inter", "assoc"])
    if len(remain) > 0:
        raise ValueError(f"Unrecognized test parameters: {remain}.")

    with session_block("QTL analysis", disable=not verbose):

        with session_line("Normalising input... ", disable=not verbose):

            data = conform_dataset(y, M, G=G)

        y = data["y"]
        M = data["M"]
        G = data["G"]
        K = data["K"]

        assert_finite(y, M, K)

        E = _asarray(E, "env", ["sample", "env"])

        if "assoc" in tests:
            slmm = StructLMM(y, M, E)
            slmm.fit(verbose=False)

        _pvi = zeros(G.shape[1])
        _pva = zeros(G.shape[1])
        for snp in range(G.shape[1]):
            x = asarray(G[:, [snp]], float)

            if "inter" in tests:
                M1 = hstack((M, x))
                slmi = StructLMM(y, M1, E)
                slmi.fit(verbose=False)
                _pvi[snp] = slmi.score_2dof_inter(x)

            if "assoc" in tests:
                _pva[snp] = slmm.score_2dof_assoc(x)

    return _pvi, _pva
