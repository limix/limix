# from __future__ import division
#
# from collections import OrderedDict
#
# import statsmodels.api as sm
# from numpy import asarray as npy_asarray
#
# from limix.util import asarray
#
# from .model import QTLModel_LM
# from .util import assure_named_covariates, named_covariates_to_array
#
#
# def qtl_test_lm(G, y, M=None, test='lrt', verbose=True):
#     """
#     Wrapper function for univariate single-variant association testing
#     using a linear model.
#
#     Args:
#         snps (ndarray):
#             (`N`, `S`) ndarray of `S` SNPs for `N` individuals.
#         pheno (ndarray):
#             (`N`, `P`) ndarray of `P` phenotype sfor `N` individuals.
#             If phenotypes have missing values, then the subset of
#             individuals used for each phenotype column will be subsetted.
#         covs (ndarray, optional):
#             (`N`, `D`) ndarray of `D` covariates for `N` individuals.
#             By default, ``covs`` is a (`N`, `1`) array of ones.
#         test ({'lrt', 'f'}, optional):
#             test statistic.
#             'lrt' for likelihood ratio test (default) or 'f' for F-test.
#         verbose (bool, optional):
#             if True, details such as runtime as displayed.
#
#     Returns:
#         :class:`limix.qtl.LMM`: LIMIX LMM object
#
#     Examples
#     --------
#
#         .. doctest::
#
#             >>> from numpy.random import RandomState
#             >>> from limix.qtl import qtl_test_lm
#             >>> from numpy import set_printoptions
#             >>> set_printoptions(4)
#             >>> random = RandomState(1)
#             >>>
#             >>> N = 100
#             >>> S = 1000
#             >>>
#             >>> snps = (random.rand(N, S) < 0.2).astype(float)
#             >>> pheno = random.randn(N, 1)
#             >>>
#             >>> lm = qtl_test_lm(snps, pheno)
#             >>> print(lm.variant_pvalues[:4])
#             [ 0.8796  0.5065  0.5666  0.6016]
#     """
#     if verbose:
#         analysis_name = "Quantitative trait locus analysis"
#         print("*** %s using LM ***" % (analysis_name, ))
#
#     G = asarray(G)
#
#     M = assure_named_covariates(M, G.shape[0])
#
#     y = npy_asarray(y, float)
#
#     lm = sm.OLS(y, named_covariates_to_array(M))
#
#     res = lm.fit()
#     null_lml = lm.loglike()
#
#     beta = res.params
#     null_covariate_effsizes = []
#
#     keys = list(M.keys())
#     for i in range(len(M)):
#         null_covariate_effsizes.append((keys[i], beta[i]))
#     null_covariate_effsizes = OrderedDict(null_covariate_effsizes)
#
#     # lm.fit(named_covariates_to_array(M), y)
#
#     # lm = qtl_test_lmm(
#     #     snps=snps, pheno=pheno, K=None, covs=covs, test=test, verbose=verbose)
#     # return lm
#
#     if verbose:
#         print(model)
#
#     return model
