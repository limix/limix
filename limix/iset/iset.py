import pandas as pd
import scipy as sp

from .mvset import MvSetTest
from .mvsetfull import MvSetTestFull
from .mvsetinc import MvSetTestInc


def fit_iSet(
    Y=None, Xr=None, F=None, Rg=None, Ug=None, Sg=None, Ie=None, n_nulls=10, factr=1e7
):
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
            the null models has iid normal residuals.
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
            (`N`, `1`) boolean ndarray indicator for analysis of
            stratified designs.
            More specifically ``Ie`` specifies in which context each
            individuals has been phenotyped.
            Needs to be specified for analysis of stratified designs.
        n_nulls (ndarray, optional):
            number of parametric bootstrap. This parameter determines
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

    Examples
    --------

        This example shows how to fit iSet when considering complete designs and
        modelling population structure/relatedness by introducing the top principle
        components of the genetic relatedness matrix (``pc_rrm``) as fixed effects.

        .. doctest::

            >>> from numpy.random import RandomState
            >>> from limix.iset import fit_iSet
            >>> from numpy import ones, concatenate
            >>> import scipy as sp
            >>>
            >>> random = RandomState(1)
            >>> sp.random.seed(0)
            >>>
            >>> N = 100
            >>> C = 2
            >>> S = 4
            >>>
            >>> snps = (random.rand(N, S) < 0.2).astype(float)
            >>> pheno = random.randn(N, C)
            >>> mean = ones((N, 1))
            >>> pc_rrm = random.randn(N, 4)
            >>> covs = concatenate([mean, pc_rrm], 1)
            >>>
            >>> df, df0 = fit_iSet(Y=pheno, Xr=snps, F=covs, n_nulls=2)
            >>>
            >>> print(df.round(3).T)  # doctest: +SKIP
                                       0
            mtSet LLR              0.166
            iSet LLR               0.137
            iSet-het LLR          -0.000
            Persistent Var         0.005
            Rescaling-GxC Var      0.005
            Heterogeneity-GxC var  0.000

        This example shows how to fit iSet when considering complete designs
        and modelling population structure/relatedness using the full
        genetic relatedness matrix.

        .. doctest::

            >>> from numpy import dot, eye
            >>> random = RandomState(1)
            >>> sp.random.seed(0)
            >>>
            >>> W = random.randn(N, 10)
            >>> kinship = dot(W, W.T) / float(10)
            >>> kinship+= 1e-4 * eye(N)
            >>>
            >>> df, df0 = fit_iSet(Y=pheno, Xr=snps, Rg=kinship, n_nulls=2)
            >>>
            >>> print(df.round(3).T)  # doctest: +SKIP
                                       0
            mtSet LLR              0.154
            iSet LLR               1.098
            iSet-het LLR           1.014
            Persistent Var         0.005
            Rescaling-GxC Var      0.005
            Heterogeneity-GxC var  0.000

        This example shows how to fit iSet when considering stratified designs
        and modelling population structure/relatedness by introducing
        the top principle components of the genetic relatedness matrix
        (``pc_rrm``) as fixed effects.
        iSet does not support models with full genetic relatedness matrix
        for stratified designs.

        .. doctest::

            >>> random = RandomState(1)
            >>> sp.random.seed(0)
            >>>
            >>> pheno = random.randn(N, 1)
            >>> Ie = random.randn(N)<0.
            >>>
            >>> df, df0 = fit_iSet(Y=pheno, Xr=snps, F=covs, Ie=Ie, n_nulls=2)
            >>>
            >>> print(df.round(3).T)  # doctest: +SKIP
                                       0
            mtSet LLR              1.177
            iSet LLR               0.648
            iSet-het LLR           0.000
            Persistent Var         0.064
            Rescaling-GxC Var      0.006
            Heterogeneity-GxC var -0.000

        For more info and examples see the `iSet tutorial`_.

        .. _iSet tutorial: https://github.com/limix/limix-tutorials/tree/master/iSet

    """
    # data
    noneNone = Sg is not None and Ug is not None
    bgRE = Rg is not None or noneNone
    # fixed effect
    msg = "The current implementation of the full rank iSet"
    msg += " does not support covariates."
    msg += " We reccommend to regress out covariates and"
    msg += " subsequently quantile normalize the phenotypes"
    msg += " to a normal distribution prior to use mtSet/iSet."
    msg += " This can be done within the LIMIX framework using"
    msg += " the methods limix.util.preprocess.regressOut and"
    msg += " limix.util.preprocess.gaussianize"
    assert not (F is not None and bgRE), msg
    # strat
    strat = Ie is not None
    msg = "iSet for interaction analysis of stratified populations "
    msg += "using contextual variables does not support random effect "
    msg += "correction for confounding. "
    msg += "Please use the fixed effects to correct for confounding. "
    assert not (strat and bgRE), msg

    # check inputs for strat design
    if strat:
        assert Y.shape[1] == 1, "Y has not the right shape"
        valid_env_ids = (sp.unique(Ie) == sp.array([0, 1])).all()
        assert valid_env_ids, "The provided Ie is not valid"

    # phenotype is centered for 2 random effect model
    if bgRE:
        Yc = Y - Y.mean(0)

    # define mtSet
    if bgRE:
        mvset = MvSetTestFull(Y=Yc, Xr=Xr, Rg=Rg, Ug=Ug, Sg=Sg, factr=factr)
    elif strat:
        mvset = MvSetTestInc(Y=Y, Xr=Xr, F=F, Ie=Ie, factr=factr)
    else:
        mvset = MvSetTest(Y=Y, Xr=Xr, F=F, factr=factr)

    RV = {}
    RV["mtSet LLR"] = mvset.assoc()
    RV["iSet LLR"] = mvset.gxe()
    RV["iSet-het LLR"] = mvset.gxehet()
    RV["Persistent Var"] = mvset.info["full"]["var_r"][0]
    RV["Rescaling-GxC Var"] = mvset.info["full"]["var_r"][1]
    RV["Heterogeneity-GxC var"] = mvset.info["full"]["var_r"][2]
    df = pd.DataFrame(RV)

    RV0 = {}
    RV0["mtSet LLR0"] = mvset.assoc_null(n_nulls=n_nulls)
    RV0["iSet LLR0"] = mvset.gxe_null(n_nulls=n_nulls)
    RV0["iSet-het LLR0"] = mvset.gxehet_null(n_nulls=n_nulls)
    df0 = pd.DataFrame(RV0)

    return df, df0
