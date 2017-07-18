from __future__ import division

from collections import OrderedDict

from numpy import asarray as npy_asarray
from numpy import diag, eye
from numpy_sugar.linalg import economic_qs

from glimix_core.glmm import GLMM
from glimix_core.lmm import LMM
from limix.stats.kinship import gower_norm
from limix.util import Timer, asarray

from .model import QTLModel
from .util import assure_named_covariates, named_covariates_to_array


def scan(G, y, lik, K=None, M=None, verbose=True):
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
        Set to ``None`` for a (generalised) linear model without random effect.
        Defaults to ``None``.
    M : (array_like, optional)
        `N` individuals by `D` covariates.
        By default, ``M`` is a (`N`, `1`) array of ones.
    verbose : (bool, optional)
        if ``True``, details such as runtime are displayed.

    Returns
    -------
    :class:`limix.qtl.model.QTLModel`
        QTL representation.

    Examples
    --------
    .. doctest::

        >>> from numpy.random import RandomState
        >>> from numpy import dot
        >>> from limix.qtl import scan
        >>> random = RandomState(1)
        >>>
        >>> N = 100
        >>> S = 1000
        >>>
        >>> snps = (random.rand(N, S) < 0.2).astype(float)
        >>> pheno = random.randn(N)
        >>> W = random.randn(N, 10)
        >>> kinship = dot(W, W.T) / float(10)
        >>>
        >>> model = scan(snps, pheno, 'normal', kinship, verbose=False)
        >>> print(model.variant_pvalues[:4])
        [ 1.      1.      0.6377  1.    ]

    .. doctest::

        >>> from numpy import dot, exp, sqrt
        >>> from numpy.random import RandomState
        >>> from limix.qtl import scan
        >>>
        >>> random = RandomState(0)
        >>>
        >>> G = random.randn(250, 500) / sqrt(500)
        >>> beta = 0.01 * random.randn(500)
        >>>
        >>> z = dot(G, beta) + 0.1 * random.randn(250)
        >>> z += dot(G[:, 0], 1) # causal SNP
        >>>
        >>> y = random.poisson(exp(z))
        >>>
        >>> candidates = G[:, :5]
        >>> K = dot(G[:, 5:], G[:, 5:].T)
        >>> model = scan(candidates, y, 'poisson', K, verbose=False)
        >>>
        >>> print(model.variant_pvalues)
        [ 0.0694  0.3336  0.5899  0.7387  0.7796]
        >>> print(model.variant_effsizes)
        [ 2.4732 -1.2588 -0.7068 -0.4772  0.3752]
        >>> print(model.variant_effsizes_se)
        [ 1.362   1.3018  1.3112  1.4309  1.3405]
    """

    if verbose:
        lik_name = lik.lower()
        lik_name = lik_name[0].upper() + lik_name[1:]
        analysis_name = "Quantitative trait locus analysis"
        print("*** %s using %s-GLMM ***" % (analysis_name, lik_name))

    G = asarray(G)
    if K is None:
        K = eye(G.shape[0])
        fix_delta = True
    else:
        fix_delta = False

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
        if fix_delta:
            lmm.delta = 1
            lmm.fix('delta')

        lmm.learn(verbose=verbose)

        null_lml = lmm.lml()

        beta = lmm.beta

        null_covariate_effsizes = []

        keys = list(M.keys())
        for i in range(len(M)):
            null_covariate_effsizes.append((keys[i], beta[i]))
        null_covariate_effsizes = OrderedDict(null_covariate_effsizes)

        flmm = lmm.get_fast_scanner()
        alt_lmls, effsizes = flmm.fast_scan(G, verbose=verbose)

        model = QTLModel(null_lml, alt_lmls, effsizes, null_covariate_effsizes)
    else:
        glmm = GLMM(y, lik, named_covariates_to_array(M), QS)
        if fix_delta:
            glmm.delta = 1
            glmm.fix('delta')
        glmm.feed().maximize(verbose=verbose)

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
        tR = diag(var - var.min() + 1e-4)
        if K is not None:
            tR += s2_g * K

        lmm = LMM(mu, X=named_covariates_to_array(M), QS=economic_qs(tR))
        lmm.learn(verbose=verbose)
        null_lml = lmm.lml()
        flmm = lmm.get_fast_scanner()
        alt_lmls, effsizes = flmm.fast_scan(G, verbose=verbose)
        model = QTLModel(null_lml, alt_lmls, effsizes, null_covariate_effsizes)

    if verbose:
        print(model)

    return model
