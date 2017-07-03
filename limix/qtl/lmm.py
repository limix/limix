from __future__ import division

from collections import OrderedDict

from glimix_core.lmm import LMM
from numpy import asarray as npy_asarray
from numpy import diag
from numpy_sugar.linalg import economic_qs

from limix.util import Timer, asarray

from .model import QTLModel_LMM
from .util import assure_named_covariates, named_covariates_to_array


def qtl_test_lmm(G, y, K, M=None, test='lrt', verbose=True):
    r"""Single-variant association testing via linear mixed models.

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
    test ({'lrt', 'f'}, optional):
        Test statistic: ``"lrt"`` for likelihood ratio test (default) or
        ``"f"`` for F-test.
    verbose : (bool, optional)
        if ``True``, details such as runtime are displayed.

    Returns
    -------
        :class:`limix.qtl.model.QTLModel_LMM`: QTL representation.

    Examples

    .. doctest::

        >>> from numpy.random import RandomState
        >>> from numpy import dot
        >>> from limix.qtl import qtl_test_lmm
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
        >>> model = qtl_test_lmm(snps, pheno, kinship, verbose=False)
        >>> print(model.variant_pvalues[:4])
        [ 1.      1.      0.6377  1.    ]
    """
    if verbose:
        analysis_name = "Quantitative trait locus analysis"
        print("*** %s using LMM ***" % (analysis_name, ))

    G = asarray(G)

    M = assure_named_covariates(M, G.shape[0])

    K = asarray(K)

    y = npy_asarray(y, float)

    desc = "Eigen decomposition of the covariance matrix..."
    with Timer(desc=desc, disable=not verbose):
        QS = economic_qs(K)

    lmm = LMM(y, X=named_covariates_to_array(M), QS=QS)
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

    model = QTLModel_LMM(null_lml, alt_lmls, effsizes, null_covariate_effsizes)

    if verbose:
        print(model)

    return model


# def forward_lmm(snps,
#                 pheno,
#                 K=None,
#                 covs=None,
#                 qvalues=False,
#                 threshold=5e-8,
#                 maxiter=2,
#                 test='lrt',
#                 verbose=None,
#                 **kw_args):
#     r"""
#     Wrapper function for univariate single-variant test with forward selection
#
#     Args:
#         snps (ndarray):
#             (`N`, `S`) ndarray of `S` SNPs for `N` individuals.
#         pheno (ndarray):
#             (`N`, `P`) ndarray of `P` phenotype sfor `N` individuals.
#             If phenotypes have missing values, then the subset of
#             individuals used for each phenotype column will be subsetted.
#         K (ndarray, optional):
#             (`N`, `N`) ndarray of LMM-covariance/kinship coefficients (optional)
#             If not provided, then standard linear regression is considered.
#         covs (ndarray, optional):
#             (`N`, `D`) ndarray of `D` covariates for `N` individuals.
#             By default, ``covs`` is a (`N`, `1`) array of ones.
#         qvalues (bool, optional):
#             if True, ``threshold`` is set to Storey qvalues to control FDR.
#             (default is False)
#         threshold (float, optional):
#             P-value thrashold for inclusion in forward selection (default 5e-8).
#             If ``qvalues=True``, the threshold is set to to Storey qvalues.
#         maxiter (int, optional):
#             maximum number of interaction scans. First scan is
#             without inclusion, so maxiter-1 inclusions can be performed.
#             (default 2)
#         test ({'lrt', 'f'}, optional):
#             test statistic.
#             'lrt' for likelihood ratio test (default) or 'f' for F-test.
#         verbose (bool, optional):
#             print verbose output. (default is False)
#
#     Returns:
#         (tuple): tuple containing:
#             - **lmm** (*:class:`limix.qtl.LMM`*): LIMIX LMM object
#             - **res** (*dict*): {**iadded**, **pvadded**, **pvall**}.
#               **iadded** is an ndarray of indices of SNPs included in order
#               of inclusion, **pvadded** is an ndarray of Pvalues obtained by
#               the included SNPs in iteration and **pvall** is a  (`Nadded`, `S`)
#               ndarray of Pvalues for all iterations
#
#     Examples
#     --------
#
#         .. doctest::
#
#             >>> from numpy.random import RandomState
#             >>> from numpy import dot, eye, ones
#             >>> from limix.qtl import forward_lmm
#             >>> random = RandomState(1)
#             >>>
#             >>> N = 100
#             >>> S = 1000
#             >>>
#             >>> snps = (random.rand(N, S) < 0.2).astype(float)
#             >>> pheno = random.randn(N, 1)
#             >>> W = random.randn(N, 10)
#             >>> kinship = dot(W, W.T) / float(10)
#             >>> kinship+= 1e-4 * eye(N)
#             >>>
#             >>> #forward lmm
#             >>> lmm, res = forward_lmm(snps, pheno, K=kinship, threshold=1.)
#             >>>
#             >>> print(res['pvall'].shape)
#             (2, 1000)
#             >>> print(res['pvall'][:,:4])
#             [[ 0.8571  0.4668  0.5872  0.5589]
#              [ 0.77    0.4226  0.6165  0.8727]]
#     """
#
#     import limix_legacy.deprecated
#     import limix_legacy.deprecated as dlimix_legacy
#     verbose = dlimix_legacy.getVerbose(verbose)
#
#     if K is None:
#         K = np.eye(snps.shape[0])
#     if covs is None:
#         covs = np.ones((snps.shape[0], 1))
#     # assert single trait
#     assert pheno.shape[1] == 1, 'forward_lmm only supports single phenotypes'
#
#     lm = qtl_test_lmm(snps, pheno.ravel(), K=K, M=covs, test=test, **kw_args)
#     pvall = []
#     pv = lm.variant_pvalues.ravel()
#     # hack to avoid issues with degenerate pv
#     pv[sp.isnan(pv)] = 1
#     pvall.append(pv)
#     imin = pv.argmin()
#     niter = 1
#     # start stuff
#     iadded = []
#     pvadded = []
#     qvadded = []
#     if qvalues:
#         assert pv.shape[
#             0] == 1, "This is untested with the fdr package. pv.shape[0]==1 failed"
#         qvall = []
#         qv = qvalues(pv)
#         qvall.append(qv)
#         score = qv.min()
#     else:
#         score = pv.min()
#     while (score < threshold) and niter < maxiter:
#         t0 = time.time()
#         iadded.append(imin)
#         pvadded.append(pv[imin])
#         if qvalues:
#             qvadded.append(qv[0, imin])
#         covs = np.concatenate((covs, snps[:, imin:(imin + 1)]), 1)
#         lm.setCovs(covs)
#         lm.process()
#         pv = lm.variant_pvalues.ravel()
#         pv[sp.isnan(pv)] = 1
#         pvall.append(pv)
#         imin = pv.argmin()
#         if qvalues:
#             qv = qvalues(pv)
#             qvall[niter:niter + 1, :] = qv
#             score = qv.min()
#         else:
#             score = pv.min()
#         t1 = time.time()
#         if verbose:
#             print(("finished GWAS testing in %.2f seconds" % (t1 - t0)))
#         niter = niter + 1
#     RV = {}
#     RV['iadded'] = iadded
#     RV['pvadded'] = pvadded
#     RV['pvall'] = np.array(pvall)
#     if qvalues:
#         RV['qvall'] = np.array(qvall)
#         RV['qvadded'] = qvadded
#     return lm, RV
