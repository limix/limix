from __future__ import division

from collections import OrderedDict

from glimix_core.glmm import GLMM
from glimix_core.lmm import LMM
from numpy import asarray as npy_asarray
from numpy import diag
from numpy_sugar.linalg import economic_qs

from limix.util import Timer, asarray

from .model import QTLModel, QTLModel_GLMM
from .util import assure_named_covariates, named_covariates_to_array


def qtl_test_glmm(G, y, lik, K, M=None, verbose=True):
    r"""Single-variant association testing via generalised linear mixed models.

    Parameters
    ----------
    G : array_like
        `N` individuals by `S` SNPs.
    y : (tuple, array_like)
        Either a tuple of two arrays of `N` individuals each (Binomial
        phenotypes) or an array of `N` individuals (Poisson or Bernoulli
        phenotypes). It does not support missing values yet.
    lik : {'bernoulli', 'binomial', 'poisson'}
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
        :class:`limix.qtl.model.QTLModel_GLMM`: QTL representation.

    Examples
    --------
    .. doctest::

        >>> from numpy import dot, exp, sqrt
        >>> from numpy.random import RandomState
        >>> from limix.qtl import qtl_test_glmm
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
        >>> model = qtl_test_glmm(candidates, y, 'poisson', K, verbose=False)
        >>>
        >>> print(model.variant_pvalues)
        [ 0.069   0.3346  0.5904  0.7392  0.779 ]
        >>> print(model.variant_effsizes)
        [ 2.4771 -1.2564 -0.706  -0.4766  0.3763]
        >>> print(model.variant_effsizes_se)
        [ 1.3624  1.3022  1.3117  1.4314  1.341 ]
    """

    if verbose:
        lik_name = lik.lower()
        lik_name = lik_name[0].upper() + lik_name[1:]
        analysis_name = "Quantitative trait locus analysis"
        print("*** %s using %s-GLMM ***" % (analysis_name, lik_name))

    G = asarray(G)

    M = assure_named_covariates(M, G.shape[0])

    K = asarray(K)

    if isinstance(y, (tuple, list)):
        y = tuple([npy_asarray(p, float) for p in y])
    else:
        y = npy_asarray(y, float)

    desc = "Eigen decomposition of the covariance matrix..."
    with Timer(desc=desc, disable=not verbose):
        QS = economic_qs(K)

    glmm = GLMM(y, lik, named_covariates_to_array(M), QS)
    glmm.feed().maximize(progress=verbose)

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
