from collections import OrderedDict

from limix._display import session_line

from .._data import conform_dataset
from .._display import session_block


def sscan(G, y, E, M=None, idx=None, tests=None, verbose=True):
    r"""
    Structured linear mixed model that accounts for genotype-environment interactions.

    StructLMM [MC18]_ extends the conventional linear mixed model by including an
    additional per-individual effect term that accounts for genotype-environment
    interaction, which can be represented as an n×1 vector, 𝛃₁.
    The model can be cast as

    .. math::

        𝐲 = 𝙼𝛂 + 𝐠𝛃₀ + 𝐠⊙𝛃₁ + 𝐞 + 𝛆, ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
        \text{where}~~ 𝛃₁∼𝓝(𝟎, 𝓋₀Σ),~~ 𝐞∼𝓝(𝟎, 𝓋₁Σ),~~\text{and}~~ 𝛆∼𝓝(𝟎, 𝓋₂𝙸).

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
    .. [MC18] Moore, R., Casale, F. P., Bonder, M. J., Horta, D., Franke, L., Barroso, I., & Stegle, O. (2018). A linear mixed-model approach to study multivariate gene–environment interactions (p. 1). Nature Publishing Group.
    """
    if tests is None:
        tests = set(["inter"])
    else:
        tests = set(tests)

    remain = tests - set(["inter", "assoc"])
    if len(remain) > 0:
        raise ValueError(f"Unrecognized test parameters: {remain}.")
    # from struct_lmm import StructLMM
    # from numpy import zeros, hstack, asarray
    # from pandas import DataFrame

    # rhos = [0.0, 0.1 ** 2, 0.2 ** 2, 0.3 ** 2, 0.4 ** 2, 0.5 ** 2, 0.5, 1.0]

    # with session_block("struct-lmm analysis", disable=not verbose):

    #     with session_line("Normalising input... ", disable=not verbose):
    #         data = conform_dataset(y, M, G=G, K=None)

    #     y = data["y"]
    #     M = data["M"]
    #     G = data["G"]

    #     if tests is None:
    #         tests = ["inter"]

    #     if "inter" in tests:
    #         slmi = StructLMM(asarray(y, float), E, W=E, rho_list=[0])

    #     if "assoc" in tests:
    #         slmm = StructLMM(asarray(y, float), E, W=E, rho_list=rhos)
    #         slmm.fit_null(F=asarray(M, float), verbose=False)

    #     _pvi = zeros(G.shape[1])
    #     _pva = zeros(G.shape[1])
    #     for snp in range(G.shape[1]):
    #         x = asarray(G[:, [snp]], float)

    #         if "inter" in tests:
    #             # interaction test
    #             M1 = hstack((M, x))
    #             slmi.fit_null(F=M1, verbose=False)
    #             _pvi[snp] = slmi.score_2_dof(x)

    #         if "assoc" in tests:
    #             # association test
    #             _pva[snp] = slmm.score_2_dof(x)

    # data = OrderedDict()
    # data["pvi"] = _pvi
    # data["pva"] = _pva
    # return DataFrame(data)
