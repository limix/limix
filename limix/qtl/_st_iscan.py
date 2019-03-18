import sys
from limix._display import session_line

from .._data import conform_dataset
from .._data import asarray as _asarray
from .._display import session_block
from ._result import ST_IScanResultFactory
from .._bits.xarray import is_dataarray
from collections import OrderedDict


def _repeat(a, repeats, axis=None):
    try:
        import dask.array as da

        return da.repeat(a, repeats, axis)
    except (TypeError, AttributeError):
        import numpy as np

        return np.repeat(a, repeats, axis)


def _cartesian_col(a, b):
    from dask.array import tile

    na = a.shape[1]
    nb = b.shape[1]

    a = tile(a, nb)
    b = _repeat(b, na, axis=1)

    if is_dataarray(a):
        a = a.data

    if is_dataarray(b):
        b = b.data

    return a * b


def _2d_sel(idx):
    from collections.abc import Iterable

    if not isinstance(idx, (slice, Iterable)):
        return [idx]

    return idx


def st_iscan(G, y, E1, idx=None, K=None, M=None, E0=None, verbose=True):
    from glimix_core.lmm import LMM as LMM
    from numpy import ones, asarray
    from numpy_sugar.linalg import economic_qs
    from xarray import concat
    from numpy_sugar import is_all_finite

    if E0 is None:
        E0 = ones([G.shape[0], 1])

    E0 = _asarray(E0, "inter0", ["sample", "inter"])
    E1 = _asarray(E1, "inter1", ["sample", "inter"])
    E01 = concat([E0, E1], dim="inter")

    data = conform_dataset(y, M, G=G, K=K)
    y = data["y"]
    M = data["M"]
    G = data["G"]
    K = data["K"]

    if idx is None:
        idx = range(G.shape[1])

    if not is_all_finite(data["y"]):
        raise ValueError("Outcome must have finite values only.")

    if not is_all_finite(data["M"]):
        raise ValueError("Covariates must have finite values only.")

    if data["K"] is not None:
        if not is_all_finite(data["K"]):
            raise ValueError("Covariate matrix must have finite values only.")
        QS = economic_qs(data["K"])
    else:
        QS = None

    lmm = LMM(y.values, M.values, QS)
    lmm.fit(verbose=verbose)
    sys.stdout.flush()

    r = ST_IScanResultFactory(M.covariate, G.candidate, E1.inter, E0.inter)
    r.set_null(lmm.lml(), lmm.beta, lmm.v1, lmm.v0)
    flmm = lmm.get_fast_scanner()
    for i in idx:

        i = _2d_sel(i)
        g = G[:, i]

        g0 = _cartesian_col(E0, g)
        g1 = _cartesian_col(E01, g)

        lml0, effs0 = flmm.scan([asarray(g0, float)], verbose=False)
        lml1, effs1 = flmm.scan([asarray(g1, float)], verbose=False)

        r.add_test(i, effs0[0], lml0[0], effs1[0], lml1[0])

    r = r.create()
    if verbose:
        print(r)

    return r


def _st_iscan(G, y, K=None, M=None, E0=None, E1=None, W_R=None, verbose=True):
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
    from numpy import asarray
    from glimix_core.lmm import LMM as LMM2
    from limix_lmm.lmm import LMM
    from limix_lmm.lmm_core import LMMCore
    from limix_core.gp import GP2KronSum, GP2KronSumLR
    from limix_core.covar import FreeFormCov
    from scipy.linalg import eigh
    from numpy import ones, var, concatenate, asarray
    from numpy_sugar.linalg import economic_qs

    lmm0 = None
    lmm2 = None
    scan = None
    scan0 = None

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
            # elif W_R is not None:
            #     if verbose:
            #         print("Model: low-rank lmm")
            #     gp = GP2KronSumLR(
            #         Y=asarray(y), Cn=FreeFormCov(1), G=W_R, F=asarray(M), A=ones((1, 1))
            #     )
            #     gp.covar.Cr.setCovariance(var(asarray(y)) * ones((1, 1)))
            #     gp.covar.Cn.setCovariance(var(asarray(y)) * ones((1, 1)))
            #     gp.optimize(verbose=verbose)
            #     Kiy_fun = gp.covar.solve
            #     QS = economic_qs(K)
            #     lmm2 = LMM2(asarray(y), asarray(M), QS)
            #     lmm2.fit(verbose=verbose)

            # case 3: full-rank linear model
            else:
                if verbose:
                    print("Model: lmm")
                # if eigh_R is None:
                eigh_R = eigh(K)
                S_R, U_R = eigh_R
                add_jitter(S_R)
                gp = GP2KronSum(
                    Y=asarray(y),
                    Cg=FreeFormCov(1),
                    Cn=FreeFormCov(1),
                    S_R=S_R,
                    U_R=U_R,
                    F=asarray(M),
                    A=ones((1, 1)),
                )
                QS = economic_qs(K)
                lmm2 = LMM2(asarray(y), asarray(M), QS)
                gp.covar.Cr.setCovariance(0.5 * var(asarray(y)) * ones((1, 1)))
                gp.covar.Cn.setCovariance(0.5 * var(asarray(y)) * ones((1, 1)))
                gp.optimize(verbose=verbose)
                lmm2.fit(verbose=verbose)
                Kiy_fun = gp.covar.solve

            if E1 is None:
                lmm = LMM(asarray(y), asarray(M), Kiy_fun)
                if lmm2 is not None:
                    scan = lmm2.get_fast_scanner()
                E1 = None
                E0 = None
            else:
                lmm = LMMCore(asarray(y), asarray(M), Kiy_fun)
                if lmm2 is not None:
                    scan = lmm2.get_fast_scanner()
                if E0 is None:
                    E0 = ones([y.shape[0], 1])
                if (E0 == 1).sum():
                    lmm0 = LMM(asarray(y), asarray(M), Kiy_fun)
                    if lmm2 is not None:
                        scan0 = lmm2.get_fast_scanner()
                else:
                    lmm0 = LMMCore(asarray(y), asarray(M), Kiy_fun)
                    if lmm2 is not None:
                        scan0 = lmm2.get_fast_scanner()
                E1 = concatenate([E0, E1], 1)

    return _process(scan, scan0, lmm, lmm0, asarray(G), E0, E1)


def _process(scan, scan0, lmm, lmm0, snps, E0, E1):
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
    from numpy import newaxis, concatenate

    if E1 is None:

        if scan is not None:
            lmls, effsizes = scan.fast_scan(snps)

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

            if scan is not None:
                lmls, effsizes = scan.fast_scan(snps)

        else:

            lmm0.process(snps, E0)

            # lmls0 = []
            # effsizes0 = []
            # if scan0 is not None:
            #     for i in range(E0.shape[1]):
            #         e0 = E0[:, i]
            #         lmlsi, effsizesi = scan.fast_scan(e0[:, newaxis] * snps)
            #         lmls0.append(lmlsi)
            #         effsizes0.append(effsizesi)

            #     lmls = concatenate(lmls)
            #     effsizes = concatenate(effsizes)

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
#     assert_likelihood(lik_name)

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
