def qtl_test_lm(snps, pheno, covs=None, test='lrt', verbose=None):
    """
    Wrapper function for univariate single-variant association testing
    using a linear model.

    Args:
        snps (ndarray):
            (`N`, `S`) ndarray of `S` SNPs for `N` individuals.
        pheno (ndarray):
            (`N`, `P`) ndarray of `P` phenotype sfor `N` individuals.
            If phenotypes have missing values, then the subset of
            individuals used for each phenotype column will be subsetted.
        covs (ndarray, optional):
            (`N`, `D`) ndarray of `D` covariates for `N` individuals.
            By default, ``covs`` is a (`N`, `1`) array of ones.
        test ({'lrt', 'f'}, optional):
            test statistic.
            'lrt' for likelihood ratio test (default) or 'f' for F-test.
        verbose (bool, optional):
            if True, details such as runtime as displayed.

    Returns:
        :class:`limix.qtl.LMM`: LIMIX LMM object

    Examples
    --------

        .. doctest::

            >>> from numpy.random import RandomState
            >>> from limix.qtl import qtl_test_lm
            >>> from numpy import set_printoptions
            >>> set_printoptions(4)
            >>> random = RandomState(1)
            >>>
            >>> N = 100
            >>> S = 1000
            >>>
            >>> snps = (random.rand(N, S) < 0.2).astype(float)
            >>> pheno = random.randn(N, 1)
            >>>
            >>> lm = qtl_test_lm(snps, pheno)
            >>> print(lm.variant_pvalues[:4])
            [ 0.8796  0.5065  0.5666  0.6016]
    """
    lm = qtl_test_lmm(
        snps=snps, pheno=pheno, K=None, covs=covs, test=test, verbose=verbose)
    return lm
