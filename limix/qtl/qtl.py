from __future__ import division

from collections import OrderedDict

from numpy import asarray as npy_asarray
from numpy import diag
from numpy_sugar.linalg import economic_qs

from glimix_core.glmm import GLMM
from glimix_core.lmm import LMM
from limix.util import Timer, asarray

from limix.stats.kinship import gower_norm
from .model import QTLModel_GLMM, QTLModel_LMM
from .util import assure_named_covariates, named_covariates_to_array


def scan(G, y, lik, K, M=None, verbose=True):
    r"""Single-variant association testing via generalised linear mixed models.

    It supports Normal, Bernoulli, Binomial, and Poisson phenotypes.
    Let :math:`N` be the sample size and :math:`S` the number of covariates.

    Parameters
    ----------
    G : array_like
        `N` individuals by `S` candidate markers.
    y : (tuple, array_like)
        Either a tuple of two arrays of `N` individuals each (Binomial
        phenotypes) or an array of `N` individuals (Normal, Poisson, or
        Bernoulli phenotypes). It does not support missing values yet.
    lik : {'normal', 'bernoulli', 'binomial', 'poisson'}
        Sample likelihood describing the residual distribution.
    K : array_like
        `N` by `N` covariance matrix (e.g., kinship coefficients).
    M : (array_like, optional)
        `N` individuals by `D` covariates.
        By default, ``M`` is a (`N`, `1`) array of ones.
    verbose : (bool, optional)
        if ``True``, details such as runtime are displayed.

    Returns
    -------
    :class:`limix.qtl.model.QTLModel`
        QTL representation.
    """

    if verbose:
        lik_name = lik.lower()
        lik_name = lik_name[0].upper() + lik_name[1:]
        analysis_name = "Quantitative trait locus analysis"
        print("*** %s using %s-GLMM ***" % (analysis_name, lik_name))

    G = asarray(G)
    K = asarray(K)
    K = gower_norm(K)
    M = assure_named_covariates(M, K.shape[0])

    if isinstance(y, (tuple, list)):
        y = tuple([npy_asarray(p, float) for p in y])
    else:
        y = npy_asarray(y, float)

    desc = "Eigen decomposition of the covariance matrix..."
    with Timer(desc=desc, disable=not verbose):
        QS = economic_qs(K)

    lik = lik.lower()

    if lik == 'normal':
        lmm = LMM(y, named_covariates_to_array(M), QS)
        lmm.learn()

        null_lml = lmm.lml()

        beta = lmm.beta

        null_covariate_effsizes = []

        keys = list(M.keys())
        for i in range(len(M)):
            null_covariate_effsizes.append((keys[i], beta[i]))
        null_covariate_effsizes = OrderedDict(null_covariate_effsizes)

        flmm = lmm.get_fast_scanner()
        alt_lmls, effsizes = flmm.fast_scan(G, verbose=verbose)

        model = QTLModel_LMM(null_lml, alt_lmls, effsizes,
                             null_covariate_effsizes)
    else:
        method = GLMM(y, lik, named_covariates_to_array(M), QS)
        method.feed().maximize(progress=verbose)

        # extract stuff from glmm
        eta = glmm._site.eta
        tau = glmm._site.tau
        scale = float(glmm.scale)
        delta = float(glmm.delta)

        beta = glmm.beta

        null_covariate_effsizes = []

        keys = list(M.keys())
        for i in range(len(M)):
            null_covariate_effsizes.append((keys[i], beta[i]))
        null_covariate_effsizes = OrderedDict(null_covariate_effsizes)

        # define useful quantities
        mu = eta / tau
        var = 1. / tau
        s2_g = scale * (1 - delta)
        tR = s2_g * K + diag(var - var.min() + 1e-4)

        lmm = LMM(mu, X=named_covariates_to_array(M), QS=economic_qs(tR))
        null_lml = lmm.lml()
        flmm = lmm.get_fast_scanner()
        alt_lmls, effsizes = flmm.fast_scan(G, verbose=verbose)
        model = QTLModel_GLMM(null_lml, alt_lmls, effsizes,
                              null_covariate_effsizes)

    if verbose:
        print(model)

    return model
