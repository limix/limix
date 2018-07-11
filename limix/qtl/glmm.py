from __future__ import division

from time import time

from numpy import asarray, diag, ones

from glimix_core.glmm import GLMMExpFam
from limix.qtl.lmm import LMM
from numpy_sugar.linalg import economic_qs


def qtl_test_glmm(
    snps,
    pheno,
    lik,
    K,
    covs=None,
    test="lrt",
    NumIntervalsDeltaAlt=100,
    searchDelta=False,
    verbose=True,
):
    """
    Wrapper function for univariate single-variant association testing
    using a generalised linear mixed model.

    Args:
        snps (array_like):
            `N` individuals by `S` SNPs.
        pheno (tuple, array_like):
            Either a tuple of two arrays of `N` individuals each (Binomial
            phenotypes) or an array of `N` individuals (Poisson or Bernoulli
            phenotypes). It does not support missing values yet.
        lik ({'bernoulli', 'binomial', 'poisson'}):
            Sample likelihood describing the residual distribution.
        K (array_like):
            `N` by `N` covariance matrix (e.g., kinship coefficients).
        covs (array_like, optional):
            `N` individuals by `D` covariates.
            By default, ``covs`` is a (`N`, `1`) array of ones.
        test ({'lrt'}, optional):
            Likelihood ratio test (default).
        NumIntervalsDeltaAlt (int, optional):
            number of steps for delta optimization on the alternative model.
            Requires ``searchDelta=True`` to have an effect.
        searchDelta (bool, optional):
            if ``True``, delta optimization on the alternative model is
            carried out. By default ``searchDelta`` is ``False``.
        verbose (bool, optional):
            if ``True``, details such as runtime are displayed.

    Returns:
        :class:`limix.qtl.LMM`: LIMIX LMM object

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
        >>> lm = qtl_test_glmm(candidates, y, 'poisson', K, verbose=False)
        >>>
        >>> print(lm.getPv())
        [[0.0694 0.3336 0.5899 0.7388 0.7796]]
    """

    snps = _asarray(snps)

    if covs is None:
        covs = ones((snps.shape[0], 1))
    else:
        covs = _asarray(covs)

    K = _asarray(K)

    if isinstance(pheno, (tuple, list)):
        y = tuple([asarray(p, float) for p in pheno])
    else:
        y = asarray(pheno, float)

    start = time()
    QS = economic_qs(K)
    glmm = GLMMExpFam(y, lik, covs, QS)
    glmm.feed().maximize(verbose=verbose)

    # extract stuff from glmm
    eta = glmm.site.eta
    tau = glmm.site.tau
    scale = float(glmm.scale)
    delta = float(glmm.delta)

    # define useful quantities
    mu = eta / tau
    var = 1. / tau
    s2_g = scale * (1 - delta)
    tR = s2_g * K + diag(var - var.min() + 1e-4)

    start = time()
    lmm = LMM(snps=snps, pheno=mu, K=tR, covs=covs, verbose=verbose)
    # if verbose:
    #     print("Elapsed time for LMM part: %.3f" % (time() - start))

    return lmm


def _asarray(X):
    import dask.array as da

    if not isinstance(X, da.Array):
        X = asarray(X, float)
    return X
