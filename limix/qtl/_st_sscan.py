from collections import OrderedDict

from limix._display import session_line

from .._data import conform_dataset
from .._display import session_block


def st_sscan(G, y, E, M=None, tests=None, verbose=True):
    """Mixed-model with genetic effect heterogeneity.

    Parameters
    ----------
    pheno : (`N`, 1) ndarray
        phenotype data
    environments : (`N`, `E`) ndarray
        environments data.
    covs : (`N`, `D`) ndarray
        covariate design matrix.
        By default, ``covs`` is a (`N`, `1`) array of ones.
    tests : list
        Which tests are performed.
        Element list values are ``'inter'`` and ``'assoc'``.
        By default, only the interaction test is considered.
    rhos : list
        for the association test, a list of ``rho`` values must be specified.
        The choice of ``rho`` affects the statistical power of the test
        (for more information see the StructLMM paper).
        By default, ``rho=[0, 0.1**2, 0.2**2, 0.3**2, 0.4**2, 0.5**2, 0.5, 1.]``
    verbose : (bool, optional):
        if True, details such as runtime as displayed.
    """
    from struct_lmm import StructLMM
    from numpy import zeros, hstack, asarray
    from pandas import DataFrame

    rhos = [0.0, 0.1 ** 2, 0.2 ** 2, 0.3 ** 2, 0.4 ** 2, 0.5 ** 2, 0.5, 1.0]

    with session_block("struct-lmm analysis", disable=not verbose):

        with session_line("Normalising input... ", disable=not verbose):
            data = conform_dataset(y, M, G=G, K=None)

        y = data["y"]
        M = data["M"]
        G = data["G"]

        if tests is None:
            tests = ["inter"]

        if "inter" in tests:
            slmi = StructLMM(asarray(y, float), E, W=E, rho_list=[0])

        if "assoc" in tests:
            slmm = StructLMM(asarray(y, float), E, W=E, rho_list=rhos)
            slmm.fit_null(F=asarray(M, float), verbose=False)

        _pvi = zeros(G.shape[1])
        _pva = zeros(G.shape[1])
        for snp in range(G.shape[1]):
            x = asarray(G[:, [snp]], float)

            if "inter" in tests:
                # interaction test
                M1 = hstack((M, x))
                slmi.fit_null(F=M1, verbose=False)
                _pvi[snp] = slmi.score_2_dof(x)

            if "assoc" in tests:
                # association test
                _pva[snp] = slmm.score_2_dof(x)

    data = OrderedDict()
    data["pvi"] = _pvi
    data["pva"] = _pva
    return DataFrame(data)
