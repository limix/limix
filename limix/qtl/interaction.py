#
# def qtl_test_interaction_lmm(snps,
#                              pheno,
#                              Inter,
#                              Inter0=None,
#                              covs=None,
#                              K=None,
#                              test='lrt'):
#     r"""
#     Wrapper function for single-variant interaction test using LMMs.
#
#     Args:
#         snps (ndarray):
#             (`N`, `S`) ndarray of `S` SNPs for `N` individuals.
#         pheno (ndarray):
#             (`N`, `P`) ndarray of `P` phenotype sfor `N` individuals.
#             If phenotypes have missing values, then the subset of
#             individuals used for each phenotype column will be subsetted.
#         Inter (ndarray, optional):
#             (`N`, `I`) ndarray of `I` interaction variables to be tested
#             for `N` individuals.
#         Inter0 (ndarray, optional):
#             (`N`, `I0`) ndarray of `I0` interaction variables to be included
#             for `N` individuals to be included both in the alternative and
#             in the null model.
#             By default `Inter0` is a (`N`, `1`) array of ones.
#         covs (ndarray, optional):
#             (`N`, `D`) ndarray of `D` covariates for `N` individuals.
#             By default, ``covs`` is a (`N`, `1`) array of ones.
#         K (ndarray, optional):
#             (`N`, `N`) ndarray of LMM-covariance/kinship coefficients (optional)
#             If not provided, then standard linear regression is considered.
#         test ({'lrt', 'f'}, optional):
#             test statistic.
#             'lrt' for likelihood ratio test (default) or 'f' for F-test.
#
#     Returns:
#         :class:`limix_legacy.deprecated.CInteractLMM()`:
#             limix CInteractLMM object
#
#     Examples
#     --------
#
#         .. doctest::
#
#             >>> from numpy.random import RandomState
#             >>> from numpy import dot, eye, ones, concatenate
#             >>> from limix.qtl import qtl_test_interaction_lmm
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
#             >>> #environments to include in null and alt models
#             >>> Inter0 = ones((N,1))
#             >>>
#             >>> #environments to include in alt model
#             >>> Inter = random.randn(N,2)
#             >>>
#             >>> #environments are also included as fixed effect covariates
#             >>> mean = ones([N, 1])
#             >>> covs = concatenate([mean, Inter], 1)
#             >>>
#             >>> #interaction test
#             >>> lmi = qtl_test_interaction_lmm(snps, pheno, Inter=Inter,
#             ...                                Inter0=Inter0)
#             >>> pvi = lmi.getPv()
#             >>>
#             >>> print(pvi.shape)
#             (1, 1000)
#             >>> print(pvi[:,:4])
#             [[ 0.8179  0.2185  0.3338  0.077 ]]
#     """
#
#     import limix_legacy.deprecated
#     import limix_legacy.deprecated as dlimix_legacy
#     N = snps.shape[0]
#     if covs is None:
#         covs = np.ones((N, 1))
#     if K is None:
#         K = np.eye(N)
#     if Inter0 is None:
#         Inter0 = np.ones([N, 1])
#     assert (pheno.shape[0] == N and K.shape[0] == N and K.shape[1] == N and
#             covs.shape[0] == N and Inter0.shape[0] == N and
#             Inter.shape[0] == N), "shapes missmatch"
#     lmi = limix_legacy.deprecated.CInteractLMM()
#     lmi.setK(K)
#     lmi.setSNPs(snps)
#     lmi.setPheno(pheno)
#     lmi.setCovs(covs)
#     lmi.setInter0(Inter0)
#     lmi.setInter(Inter)
#     if test == 'lrt':
#         lmi.setTestStatistics(lmi.TEST_LRT)
#     elif test == 'f':
#         lmi.setTestStatistics(lmi.TEST_F)
#     else:
#         print(test)
#         raise NotImplementedError("only f or lrt are implemented")
#     lmi.process()
#     return lmi
#
#
#
# # TODO: we need to fix. THis does not work as interact_GxE is not existing
# # I vote we also use **kw_args to forward parameters to interact_Gxe?
# def qtl_test_interaction_GxG(pheno,
#                              snps1,
#                              snps2=None,
#                              K=None,
#                              covs=None,
#                              test='lrt'):
#     """
#     Epistasis test between two sets of SNPs
#
#     Args:
#         pheno:  [N x 1] np.array of 1 phenotype for N individuals
#         snps1:  [N x S1] np.array of S1 SNPs for N individuals
#         snps2:  [N x S2] np.array of S2 SNPs for N individuals
#         K:      [N x N] np.array of LMM-covariance/kinship koefficients (optional)
#                         If not provided, then linear regression analysis is performed
#         covs:   [N x D] np.array of D covariates for N individuals
#         test:    'lrt' for likelihood ratio test (default) or 'f' for F-test
#
#     Returns:
#         pv:     [S2 x S1] np.array of P values for epistasis tests beten all SNPs in
#                 snps1 and snps2
#     """
#     if K is None:
#         K = np.eye(N)
#     N = snps1.shape[0]
#     if snps2 is None:
#         snps2 = snps1
#     return qtl_test_interaction_GxE_1dof(
#         snps=snps1, pheno=pheno, env=snps2, covs=covs, K=K, test=test)
#
#
# def qtl_test_interaction_GxE_1dof(snps,
#                                   pheno,
#                                   env,
#                                   K=None,
#                                   covs=None,
#                                   test='lrt',
#                                   verbose=None):
#     """
#     Univariate GxE fixed effects interaction linear mixed model test for all
#     pairs of SNPs and environmental variables.
#
#     Args:
#         snps:   [N x S] np.array of S SNPs for N individuals
#         pheno:  [N x 1] np.array of 1 phenotype for N individuals
#         env:    [N x E] np.array of E environmental variables for N individuals
#         K:      [N x N] np.array of LMM-covariance/kinship koefficients (optional)
#                         If not provided, then linear regression analysis is performed
#         covs:   [N x D] np.array of D covariates for N individuals
#         test:    'lrt' for likelihood ratio test (default) or 'f' for F-test
#         verbose: print verbose output? (False)
#
#     Returns:
#         pv:     [E x S] np.array of P values for interaction tests between all
#                 E environmental variables and all S SNPs
#     """
#     import limix_legacy.deprecated as dlimix_legacy
#     verbose = dlimix_legacy.getVerbose(verbose)
#     N = snps.shape[0]
#     if K is None:
#         K = np.eye(N)
#     if covs is None:
#         covs = np.ones((N, 1))
#     assert (env.shape[0] == N and pheno.shape[0] == N and K.shape[0] == N and
#             K.shape[1] == N and covs.shape[0] == N), "shapes missmatch"
#     Inter0 = np.ones((N, 1))
#     pv = np.zeros((env.shape[1], snps.shape[1]))
#     if verbose:
#         print(("starting %i interaction scans for %i SNPs each." %
#                (env.shape[1], snps.shape[1])))
#     t0 = time.time()
#     for i in range(env.shape[1]):
#         t0_i = time.time()
#         cov_i = np.concatenate((covs, env[:, i:(i + 1)]), 1)
#         lm_i = qtl_test_interaction_lmm(
#             snps=snps,
#             pheno=pheno,
#             covs=cov_i,
#             Inter=env[:, i:(i + 1)],
#             Inter0=Inter0,
#             test=test)
#         pv[i, :] = lm_i.getPv()[0, :]
#         t1_i = time.time()
#         if verbose:
#             print(("Finished %i out of %i interaction scans in %.2f seconds." %
#                    ((i + 1), env.shape[1], (t1_i - t0_i))))
#     t1 = time.time()
#     print((
#         "-----------------------------------------------------------\nFinished all %i interaction scans in %.2f seconds."
#         % (env.shape[1], (t1 - t0))))
#     return pv
