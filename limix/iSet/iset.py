from .mvSet import MvSetTest
from .mvSetFull import MvSetTestFull
from .mvSetInc import MvSetTestInc
import pandas as pd

def fit_iSet(Y=None, Xr=None, F=None, Rg=None, Ug=None, Sg=None, Ie=None, n_nulls=10, factr=1e7):
    """
    Fit interaction set test (iSet).

    Args:
        Y (ndarray):
            For complete design, the phenotype ndarray `Y` 
            for `N` samples and `C` contexts has shape (`N`, `C`).
            For stratified design, the phenotype ndarray `Y` has
            shape (`N`, `1`) (each individual is phenotyped in only one context -
            ``Ie`` specifies in which context each individuals has been phenotyped).
        Xr (ndarray):
            (`N`, `S`) genotype values for `N` samples and `S` variants
            (defines the set component)
        F (ndarray, optional):
            (`N`, `K`) ndarray of `K` covariates for `N` individuals.
            By default, ``F`` is a (`N`, `1`) array of ones.
        Rg (ndarray, optional):
            (`N`, `N`) ndarray of LMM-covariance/kinship coefficients.
            ``Ug`` and ``Sg`` can be provided instead of ``Rg``.
            If neither ``Rg`` nor ``Ug`` and ``Sg`` are provided,
            iid normal residuals are considered.
        Ug (ndarray, optional):
            (`N`, `N`) ndarray of eigenvectors of ``Rg``. 
            ``Ug`` and ``Sg`` can be provided instead of ``Rg``.
            If neither ``Rg`` nor ``Ug`` and ``Sg`` are provided,
            iid normal residuals are considered.
        Sg (ndarray, optional):
            (`N`, ) ndarray of eigenvalues of ``Rg``.
            ``Ug`` and ``Sg`` can be provided instead of ``Rg``.
            If neither ``Rg`` nor ``Ug`` and ``Sg`` are provided,
            iid normal residuals are considered.
        Ie (ndarry, optional):
            (`N`, `1`) binary indicator for analysis of stratified designs.
            More specifically ``Ie`` specifies in which context each
            individuals has been phenotyped. 
        n_nulls (ndarray, optional):
            number of parametric bootstrap. This parameter determine
            the minimum P value that can be estimated.
            The default value is 10.
        factr (float, optional):
            optimization paramenter that determines the accuracy of the solution
            (see scipy.optimize.fmin_l_bfgs_b for more details).

    Returns:
        (tuple): tuple containing:
            - **df** (*:class:`pandas.DataFrame`*):
              contains test statistcs of mtSet, iSet, and iSet-GxC tests and
              the variance exaplained by persistent, GxC and
              heterogeneity-GxC effects.
            - **df0** (*:class:`pandas.DataFrame`*):
              contains null test statistcs of mtSet, iSet, and iSet-GxC tests.
    """
    # TODO: add condition C==2
    # TODO: add condition that Ie and Rg do not co-occurr
    # TODO: add conditions in a way that the code does not break

    # data
    noneNone = Sg is not None and Ug is not None
    bgRE = Rg is not None or noneNone
    # fixed effect
    msg = 'The current implementation of the full rank iSet'
    msg+= ' does not support covariates.'
    msg+= ' We reccommend to regress out covariates and'
    msg+= ' subsequently quantile normalize the phenotypes'
    msg+= ' to a normal distribution prior to use mtSet/iSet.'
    msg+= ' This can be done within the LIMIX framework using'
    msg+= ' the methods limix.util.preprocess.regressOut and'
    msg+= ' limix.util.preprocess.gaussianize'
    assert not (F is not None and bgRE), msg
    # strat
    strat = Ie is not None
    msg = 'iSet for interaction analysis of stratified populations ' 
    msg+= 'using contextual variables does not support random effect '
    msg+= 'correction for confounding. '
    msg+= 'Please use the fixed effects to correct for confounding. '
    assert not (strat and bgRE), msg

    #define mtSet
    if bgRE:        mvset = MvSetTestFull(Y=Y,Xr=Xr,Rg=Rg,Ug=Ug,Sg=Sg,factr=factr)
    elif strat:     mvset = MvSetTestInc(Y=Y,Xr=Xr,F=F,Ie=Ie,factr=factr)
    else:           mvset = MvSetTest(Y=Y,Xr=Xr,F=F,factr=factr)

    RV = {}
    RV['mtSet LLR'] = mvset.assoc()
    RV['iSet LLR'] = mvset.gxe()
    RV['iSet-het LLR'] = mvset.gxehet()
    RV['Persistent Var'] = mvset.info['full']['var_r'][0]
    RV['Rescaling-GxC Var'] = mvset.info['full']['var_r'][1]
    RV['Heterogeneity-GxC var'] = mvset.info['full']['var_r'][2]
    df = pd.DataFrame(RV)

    RV0 = {}
    RV0['mtSet LLR0'] = mvset.assoc_null(n_nulls=n_nulls)
    RV0['iSet LLR0'] = mvset.gxe_null(n_nulls=n_nulls)
    RV0['iSet-het LLR0'] = mvset.gxehet_null(n_nulls=n_nulls)
    df0 = pd.DataFrame(RV0)

    return df, df0

