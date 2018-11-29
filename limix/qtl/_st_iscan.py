from limix._display import session_line

from .._data import conform_dataset
from .._display import session_block
from collections import OrderedDict


def st_iscan(G, y, K=None, M=None, E0=None, E1=None, W_R=None, verbose=True):
    r""" Single-variant association interation testing.

    Parameters
    ----------
    pheno : (`N`, 1) ndarray
        phenotype data
    covs : (`N`, `D`) ndarray
        covariate design matrix.
        By default, ``covs`` is a (`N`, `1`) array of ones.
    R : (`N`, `N`) ndarray
        LMM-covariance/genetic relatedness matrix.
        If not provided, then standard linear regression is considered.
        Alternatively, its eighenvalue decomposition can be
        provided through ``eigh_R``.
        if ``eigh_R`` is set, this parameter is ignored.
        If the LMM-covariance is low-rank, ``W_R`` can be provided
    eigh_R : tuple
        Tuple with `N` ndarray of eigenvalues of `R` and
        (`N`, `N`) ndarray of eigenvectors of ``R``.
    W_R : (`N`, `R`) ndarray
        If the LMM-covariance is low-rank, one can provide ``W_R`` such that
        ``R`` = dot(``W_R``, transpose(``W_R``)).
    inter : (`N`, `K`) ndarray
        interaction variables interacting with the snp.
        If specified, then the current tests are considered:
        (i) (inter&inter0)-by-g vs no-genotype-effect;
        (ii) inter0-by-g vs no-genotype-effect;
        (iii) (inter&inter0)-by-g vs inter0-by-g.
    inter0 : (`N`, `K0`) ndarray
        interaction variables to be included in the alt and null model.
        By default, if inter is not specified, inter0 is ignored.
        By default, if inter is specified, inter0=ones so that inter0-by-g=g,
        i.e. an additive genetic effect is considered.
    verbose : (bool, optional):
        if True, details such as runtime as displayed.
    """
    from limix_lmm.lmm import LMM
    from limix_lmm.lmm_core import LMMCore
    from limix_core.gp import GP2KronSum, GP2KronSumLR
    from limix_core.covar import FreeFormCov
    from scipy.linalg import eigh
    from numpy import ones, var, concatenate, asarray

    lmm0 = None

    with session_block("single-trait association test", disable=not verbose):

        # if covs is None:
        #     covs = ones([pheno.shape[0], 1])

        with session_line("Normalising input... ", disable=not verbose):
            data = conform_dataset(y, M, G=G, K=K)

            y = data["y"]
            M = data["M"]
            G = data["G"]
            K = data["K"]

            # case 1: linear model
            # if W_R is None and eigh_R is None and R is None:
            if K is None:
                if verbose:
                    print("Model: lm")
                gp = None
                Kiy_fun = None

            # case 2: low-rank linear model
            elif W_R is not None:
                if verbose:
                    print("Model: low-rank lmm")
                gp = GP2KronSumLR(Y=y, Cn=FreeFormCov(1), G=W_R, F=M, A=ones((1, 1)))
                gp.covar.Cr.setCovariance(var(y) * ones((1, 1)))
                gp.covar.Cn.setCovariance(var(y) * ones((1, 1)))
                gp.optimize(verbose=verbose)
                Kiy_fun = gp.covar.solve

            # case 3: full-rank linear model
            else:
                if verbose:
                    print("Model: lmm")
                # if eigh_R is None:
                eigh_R = eigh(K)
                S_R, U_R = eigh_R
                add_jitter(S_R)
                gp = GP2KronSum(
                    Y=y,
                    Cg=FreeFormCov(1),
                    Cn=FreeFormCov(1),
                    S_R=S_R,
                    U_R=U_R,
                    F=M,
                    A=ones((1, 1)),
                )
                gp.covar.Cr.setCovariance(0.5 * var(y) * ones((1, 1)))
                gp.covar.Cn.setCovariance(0.5 * var(y) * ones((1, 1)))
                gp.optimize(verbose=verbose)
                Kiy_fun = gp.covar.solve

            if E1 is None:
                lmm = LMM(y, M, Kiy_fun)
                E1 = None
                E0 = None
            else:
                lmm = LMMCore(y, M, Kiy_fun)
                if E0 is None:
                    E0 = ones([y.shape[0], 1])
                if (E0 == 1).sum():
                    lmm0 = LMM(y, M, Kiy_fun)
                else:
                    lmm0 = LMMCore(y, M, Kiy_fun)
                E1 = concatenate([E0, E1], 1)

    return _process(lmm, lmm0, asarray(G), E0, E1)


def _process(lmm, lmm0, snps, E0, E1):
    """
    Parameters
    ----------
    snps : (`N`, `S`) ndarray
        genotype data
    return_ste : bool
        if True, return eff size standard errors(default is False)
    return_lrt : bool
        if True, return llr test statistcs (default is False)

    Return
    ------
    res : pandas DataFrame
        Results as pandas dataframs
    """
    from scipy.stats import chi2
    from pandas import DataFrame

    if E1 is None:

        lmm.process(snps)
        RV = OrderedDict()
        RV["pv"] = lmm.getPv()
        RV["beta"] = lmm.getBetaSNP()
        RV["beta_ste"] = lmm.getBetaSNPste()
        RV["lrt"] = lmm.getLRT()

    else:

        lmm.process(snps, E1)
        if (E0 == 1).sum():
            lmm0.process(snps)
        else:
            lmm0.process(snps, E0)

        # compute pv
        lrt1 = lmm.getLRT()
        lrt0 = lmm0.getLRT()
        lrt = lrt1 - lrt0
        pv = chi2(E1.shape[1] - E0.shape[1]).sf(lrt)

        RV = OrderedDict()
        RV["pv1"] = lmm.getPv()
        RV["pv0"] = lmm0.getPv()
        RV["pv"] = pv
        if (E0 == 1).sum():
            RV["beta0"] = lmm0.getBetaSNP()
            RV["beta0_ste"] = lmm0.getBetaSNPste()
        RV["lrt1"] = lrt1
        RV["lrt0"] = lrt0
        RV["lrt"] = lrt

    return DataFrame(RV)


# def st_iscan(G, y, lik, inter, Ginter=None, K=None, M=None, verbose=True):
#     r"""Interaction single-variant association testing via mixed models.

#     It supports Normal (linear mixed model), Bernoulli, Binomial, and Poisson
#     residual errors, defined by ``lik``.
#     The columns of ``G`` define the candidates to be tested for association
#     with the phenotype ``y``.
#     The covariance matrix is set by ``K``.
#     If not provided, or set to ``None``, the generalised linear model
#     without random effects is assumed.
#     The covariates can be set via the parameter ``M``.
#     We recommend to always provide a column of ones in the case

#     Parameters
#     ----------
#     G : array_like
#         `n` individuals by `s` candidate markers.
#     y : tuple, array_like
#         Either a tuple of two arrays of `n` individuals each (Binomial
#         phenotypes) or an array of `n` individuals (Normal, Poisson, or
#         Bernoulli phenotypes). It does not support missing values yet.
#     lik : {'normal', 'bernoulli', 'binomial', 'poisson'}
#         Sample likelihood describing the residual distribution.
#     inter : array_like
#         `n` individuals by `i` interaction factors.
#     Ginter : array_like
#         `n` individuals by `s` candidate markers. used for interaction.
#         Defaults to ``None``, in which case G is used.
#     K : array_like, optional
#         `n` by `n` covariance matrix (e.g., kinship coefficients).
#         Set to ``None`` for a generalised linear model without random effects.
#         Defaults to ``None``.
#     M : array_like, optional
#         `n` individuals by `d` covariates.
#         By default, ``M`` is an `n` by `1` matrix of ones.
#     verbose : bool, optional
#         if ``True``, details such as runtime are displayed.

#     Returns
#     -------
#     :class:`limix.qtl.model.QTLModel`
#         Interaction QTL representation.
#     """
#     from numpy_sugar import is_all_finite
#     from numpy_sugar.linalg import economic_qs

#     lik = lik.lower()

#     if not isinstance(lik, (tuple, list)):
#         lik = (lik,)

#     lik_name = lik[0].lower()
#     check_likelihood_name(lik_name)

#     if Ginter is None:
#         Ginter = G

#     with session_block("interaction qtl analysis", disable=not verbose):

#         with session_line("Normalising input... ", disable=not verbose):
#             data = conform_dataset(
#                 y,
#                 M,
#                 G=G,
#                 K=K,
#                 X=[
#                     (inter, "interactions", "interaction"),
#                     (Ginter, "icandidates", "icandidate"),
#                 ],
#             )

#         y = data["y"]
#         M = data["M"]
#         G = data["G"]
#         K = data["K"]
#         inter = data["X0"]
#         Ginter = data["X1"]

#         if not is_all_finite(y):
#             raise ValueError("Outcome must have finite values only.")

#         if not is_all_finite(M):
#             raise ValueError("Covariates must have finite values only.")

#         if K is not None:
#             if not is_all_finite(K):
#                 raise ValueError("Covariate matrix must have finite values only.")
#             QS = economic_qs(K)
#         else:
#             QS = None

#         y = normalise_extreme_values(data["y"], lik)

#         if lik_name == "normal":
#             model = _perform_lmm(y.values, M, QS, G, inter, Ginter, verbose)
#         else:
#             model = _perform_glmm(y.values, lik, M, K, QS, G, inter, Ginter, verbose)

#         return model


# def _perform_lmm(y, M, QS, G, inter, Ginter, verbose):
#     from xarray import DataArray
#     from glimix_core.lmm import LMM
#     from tqdm import tqdm
#     from numpy import concatenate, newaxis

#     alt_lmls = []
#     effsizes = []
#     ncov_effsizes = []
#     null_lmls = []
#     interv = inter.values

#     for i in tqdm(range(Ginter.shape[1]), disable=not verbose):
#         gi = Ginter[:, i].values[:, newaxis]
#         X1 = gi * interv

#         covariates = concatenate((M.values, gi), axis=1)
#         lmm = LMM(y, covariates, QS)

#         lmm.fit(verbose=verbose)

#         null_lmls.append(lmm.lml())

#         ncov_effsizes.append(lmm.beta)

#         flmm = lmm.get_fast_scanner()
#         alt_lmls_, effsizes_ = flmm.fast_scan(X1, verbose=verbose)
#         alt_lmls.append(alt_lmls_)
#         effsizes.append(effsizes_)

#     alt_lmls = DataArray(
#         data=alt_lmls,
#         dims=["interaction", "candidate"],
#         coords={
#             "interaction": inter.coords["interaction"],
#             "candidate": G.coords["candidate"],
#         },
#     )
#     effsizes = DataArray(data=effsizes, dims=["interaction", "candidate"])

#     index = list(M.columns) + ["variant"]
#     ncov_effsizes = DataArray(data=ncov_effsizes, index=index)
#     null_lmls = DataArray(
#         null_lmls,
#         dims=["interaction"],
#         coords={"interaction": inter.coords["interaction"]},
#     )

#     return QTLModel(null_lmls, alt_lmls, effsizes, ncov_effsizes)


# def _perform_glmm(y, M, QS, G, inter, Ginter, verbose):
#     raise NotImplementedError


# Examples
# --------
# .. doctest::

#     >>> from numpy import dot, zeros, asarray
#     >>> from numpy.random import RandomState
#     >>> from numpy.testing import assert_allclose
#     >>>
#     >>> from limix.qtl import st_iscan
#     >>> from pandas import DataFrame, option_context
#     >>>
#     >>> random = RandomState(0)
#     >>> nsamples = 50
#     >>>
#     >>> X = random.randn(nsamples, 10)
#     >>> G = random.randn(nsamples, 100)
#     >>> K = dot(G, G.T)
#     >>> ntrials = random.randint(1, 100, nsamples)
#     >>> z = dot(G, random.randn(100)) / 10
#     >>>
#     >>> successes = zeros(len(ntrials), int)
#     >>> for i in range(len(ntrials)):
#     ...     for j in range(ntrials[i]):
#     ...         successes[i] += int(z[i] + 0.5 * random.randn() > 0)
#     >>>
#     >>> y = successes / asarray(ntrials, float)
#     >>>
#     >>> inter = random.randn(nsamples, 3)
#     >>>
#     >>> index = ['sample%02d' % i for i in range(X.shape[0])]
#     >>> cols = ['SNP%02d' % i for i in range(X.shape[1])]
#     >>> X = DataFrame(data=X, index=index, columns=cols)
#     >>>
#     >>> cols = ['inter%02d' % i for i in range(inter.shape[1])]
#     >>> inter = DataFrame(data=inter, index=index, columns=cols)
#     >>>
#     >>> model = st_iscan(X, y, 'normal', inter, K, verbose=False)
#     >>>
#     >>> with option_context('precision', 5):
#     ...     print(model.variant_pvalues)
#            inter00  inter01  inter02
#     SNP00  0.81180  0.63035  0.61240
#     SNP01  0.02847  0.64437  0.82671
#     SNP02  0.56817  0.72882  0.23928
#     SNP03  0.53793  0.64628  0.86144
#     SNP04  0.13858  0.39475  0.28650
#     SNP05  0.06722  0.56295  0.39859
#     SNP06  0.12739  0.62219  0.68084
#     SNP07  0.32834  0.96894  0.67628
#     SNP08  0.28341  0.29361  0.56248
#     SNP09  0.64945  0.67185  0.76600


def add_jitter(S_R):
    from numpy import maximum

    assert S_R.min() > -1e-6, "LMM-covariance is not sdp!"
    RV = S_R + maximum(1e-4 - S_R.min(), 0)
    return RV
