# Copyright(c) 2014, The LIMIX developers (Christoph Lippert, Paolo Francesco Casale, Oliver Stegle)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
qtl.py contains wrappers around C++ Limix objects to streamline common tasks in GWAS.
"""

import time

import numpy as np
import scipy as sp
import scipy.stats as st

from limix.qtl.lmm import LMM
from limix.stats import qvalues


def qtl_test_lm(snps, pheno, covs=None, test="lrt", verbose=None):
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
            >>> print(lm.getPv()[:,:4])
            [[0.8796 0.5065 0.5666 0.6016]]
    """
    lm = qtl_test_lmm(
        snps=snps, pheno=pheno, K=None, covs=covs, test=test, verbose=verbose
    )
    return lm


def qtl_test_lmm(
    snps,
    pheno,
    K=None,
    covs=None,
    test="lrt",
    NumIntervalsDelta0=100,
    NumIntervalsDeltaAlt=100,
    searchDelta=False,
    verbose=None,
):
    """
    Wrapper function for univariate single-variant association testing
    using a linear mixed model.

    Args:
        snps (ndarray):
            (`N`, `S`) ndarray of `S` SNPs for `N` individuals.
        pheno (ndarray):
            (`N`, `P`) ndarray of `P` phenotype sfor `N` individuals.
            If phenotypes have missing values, then the subset of
            individuals used for each phenotype column will be subsetted.
        K (ndarray, optional):
            (`N`, `N`) ndarray of LMM-covariance/kinship coefficients.
            If not provided, then standard linear regression is considered.
        covs (ndarray, optional):
            (`N`, `D`) ndarray of `D` covariates for `N` individuals.
            By default, ``covs`` is a (`N`, `1`) array of ones.
        test ({'lrt', 'f'}, optional):
            test statistic.
            'lrt' for likelihood ratio test (default) or 'f' for F-test.
        NumIntervalsDelta0 (int, optional):
            number of steps for delta optimization on the null model.
            By default ``NumIntervalsDelta0`` is 100.
        NumIntervalsDeltaAlt (int, optional):
            number of steps for delta optimization on the alternative model.
            Requires ``searchDelta=True`` to have an effect.
        searchDelta (bool, optional):
            if True, delta optimization on the alternative model is carried out.
            By default ``searchDelta`` is False.
        verbose (bool, optional):
            if True, details such as runtime as displayed.

    Returns:
        :class:`limix.qtl.LMM`: LIMIX LMM object

    Examples
    --------

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
            >>> pheno = random.randn(N, 1)
            >>> W = random.randn(N, 10)
            >>> kinship = dot(W, W.T) / float(10)
            >>>
            >>> lmm = qtl_test_lmm(snps, pheno, kinship)
            >>> print(lmm.getPv()[:,:4])
            [[0.8571 0.4668 0.5872 0.5589]]
    """
    lmm_ = LMM(
        snps=snps,
        pheno=pheno,
        K=K,
        covs=covs,
        test=test,
        NumIntervalsDelta0=NumIntervalsDelta0,
        NumIntervalsDeltaAlt=NumIntervalsDeltaAlt,
        searchDelta=searchDelta,
        verbose=verbose,
    )
    return lmm_


def qtl_test_lmm_kronecker(
    snps,
    phenos,
    covs=None,
    Acovs=None,
    Asnps=None,
    K1r=None,
    K1c=None,
    K2r=None,
    K2c=None,
    NumIntervalsDelta0=100,
    NumIntervalsDeltaAlt=100,
    searchDelta=False,
):
    r"""
    Wrapper function for single-variant multi-trait association testing class
    using linear mixed models.

    Args:
        snps (ndarray):
            (`N`, `S`) ndarray of `S` SNPs for `N` individuals.
        pheno (ndarray):
            (`N`, `P`) ndarray of `P` phenotype sfor `N` individuals.
            If phenotypes have missing values, then the subset of
            individuals used for each phenotype column will be subsetted.
        covs (ndarray or list of ndarrays, optional):
            sample design matrices of the different fixed effect terms.
            In the case of a single fixed effect term, ``covs`` can be specified
            as a (`N`, `D`) ndarray of `D` covariates for `N` individuals.
            In this case ``Acovs`` should also be a (`K`, `P`) ndarray,
            where `K` is the rank of the trait design and
            `P` is the number of traits.
            For multiple terms, ``covs`` is a list of ndarray, where each
            element ``covs[i]`` has shape (`N`, `Di`).
            In this case, ``Acovs`` should also be list with the same number
            of  ndarray elements.
            Each element ``Acovs[i]`` has shape (`Ki`, `P`), where `P` is the
            number of analyzed phenotypes.
            By default covs is (`N`, `1`) ndarray of ones and
            ``Acovs[i]`` is a (`P`, `P`) identity ndarray.
        Acovs (ndarray or list of ndarrays, optional):
            trait design matrices of the different fixed effect terms.
            See `covs` for more info.
        Asnps (ndarray, optional):
            (`K`, `P`) trait design of rank `K` of the snp fixed effect on the
            `P` traits. By default ``Asnps`` is (`1`, `P`) ndrray of ones.
        K1r (ndarray, optional):
            (`N`, `N`) ndarray of LMM-covariance/kinship coefficients.
            If not provided, ``K1r = numpy.dot(snps, snps.T)`` is assumed.
        K1c (ndarray, optional):
            (`P`, `P`) ndarray of LMM-covariance/kinship coefficients.
            If either `K1c` or `K2c` is not provided, K1c is estimated by
            (restricted) maximum likelihood  from the model without the snp effect.
        K2r (ndarray, optional):
            (`N`, `N`) ndarray of LMM-covariance/kinship coefficients.
            By default ``K2r`` is a (`N`, `N`) identity ndarray.
        K2c (ndarray, optional):
            (`P`, `P`) ndarray of residual covariances.
            If either `K1c` or `K2c` is not provided, K1c is estimated by
            (restricted) maximum likelihood  from the model without the snp effect.
        NumIntervalsDelta0 (int, optional):
            number of steps for delta optimization on the null model.
            By default ``NumIntervalsDelta0`` is 100.
        NumIntervalsDeltaAlt (int, optional):
            number of steps for delta optimization on the alternative model.
            Requires ``searchDelta=True`` to have an effect.
        searchDelta (bool, optional):
            if True, delta optimization on the alternative model is carried out.
            By default ``searchDelta`` is False.

    Returns:
        (tuple): tuple containing:
            - **lmm** (*:class:`limix_legacy.deprecated.CKroneckerLMM()`*):
              CKroneckerLMM limix_legacy object
            - **pv** (*ndarray*):
              P-values for all SNPs from liklelihood ratio test

    Examples
    --------

        Example of multi-trait single-variant association testing using
        a linear mixed model and specifying a common effects.

        .. doctest::

            >>> from numpy.random import RandomState
            >>> from numpy import dot, eye, ones
            >>> from limix.qtl import qtl_test_lmm_kronecker
            >>> random = RandomState(1)
            >>>
            >>> N = 100
            >>> S = 1000
            >>> P = 2
            >>>
            >>> snps = (random.rand(N, S) < 0.2).astype(float)
            >>> pheno = random.randn(N, P)
            >>> W = random.randn(N, 10)
            >>> kinship = dot(W, W.T) / float(10)
            >>> kinship+= 1e-4 * eye(N)
            >>>
            >>> #common effect test
            >>> Asnps = ones((1, P))
            >>> lmm, pv = qtl_test_lmm_kronecker(snps, pheno, K1r=kinship,
            ...                                 Asnps=Asnps)
            >>> beta = lmm.getBetaSNP()
            >>>
            >>> print(pv.shape)
            (1, 1000)
            >>> print(pv[:,:4])
            [[0.8722 0.6509 0.3593 0.6816]]
            >>>
            >>> print(beta.shape)
            (1, 1000)
            >>> print(beta[:,:4])
            [[0.0281 0.0758 0.1458 0.0701]]

        Example showing how to set an any effect test.
        For more information on effect designs
        see `limix_tutorials.`_

        .. _limix_tutorials.: http://nbviewer.jupyter.org/github/limix/limix-tutorials/blob/master/index.ipynb

        .. doctest::

            >>> #any effect test
            >>> Asnps = eye(P)
            >>> lmm, pv = qtl_test_lmm_kronecker(snps, pheno, K1r=kinship,
            ...                                 Asnps=Asnps)
            >>> beta = lmm.getBetaSNP()
            >>>
            >>> print(pv.shape)
            (1, 1000)
            >>> print(pv[:,:4])
            [[0.9838 0.8664 0.3555 0.5244]]
            >>>
            >>> print(beta.shape)
            (2, 1000)
            >>> print(beta[:,:4]) # doctest: +SKIP
            [[ 0.0413  0.0312 -0.0177 -0.098 ]
             [ 0.0139  0.1234  0.3204  0.2495]]
    """
    import limix_legacy.deprecated
    import limix_legacy.deprecated as dlimix_legacy

    # 0. checks
    N = phenos.shape[0]
    P = phenos.shape[1]

    if K1r is None:
        K1r = np.dot(snps, snps.T)
    else:
        assert K1r.shape[0] == N, "K1r: dimensions dismatch"
        assert K1r.shape[1] == N, "K1r: dimensions dismatch"

    if K2r is None:
        K2r = np.eye(N)
    else:
        assert K2r.shape[0] == N, "K2r: dimensions dismatch"
        assert K2r.shape[1] == N, "K2r: dimensions dismatch"

    covs, Acovs = _updateKronCovs(covs, Acovs, N, P)

    # Asnps can be several designs
    if Asnps is None:
        Asnps = [np.ones([1, P])]
    if not isinstance(Asnps, list):
        Asnps = [Asnps]
    assert len(Asnps) > 0, "need at least one Snp design matrix"

    # one row per column design matrix
    pv = np.zeros((len(Asnps), snps.shape[1]))

    # 1. run GP model to infer suitable covariance structure
    if K1c is None or K2c is None:
        vc = _estimateKronCovariances(
            phenos=phenos,
            K1r=K1r,
            K2r=K2r,
            K1c=K1c,
            K2c=K2c,
            covs=covs,
            Acovs=Acovs,
            verbose=False,
        )
        K1c = vc.getTraitCovar(0)
        K2c = vc.getTraitCovar(1)
    else:
        assert K1c.shape[0] == P, "K1c: dimensions dismatch"
        assert K1c.shape[1] == P, "K1c: dimensions dismatch"
        assert K2c.shape[0] == P, "K2c: dimensions dismatch"
        assert K2c.shape[1] == P, "K2c: dimensions dismatch"

    # 2. run kroneckerLMM
    lmm = limix_legacy.deprecated.CKroneckerLMM()
    lmm.setK1r(K1r)
    lmm.setK1c(K1c)
    lmm.setK2r(K2r)
    lmm.setK2c(K2c)
    lmm.setSNPs(snps)
    # add covariates
    for ic in range(len(Acovs)):
        lmm.addCovariates(covs[ic], Acovs[ic])
    lmm.setPheno(phenos)

    # delta serch on alt. model?
    if searchDelta:
        lmm.setNumIntervalsAlt(NumIntervalsDeltaAlt)
    else:
        lmm.setNumIntervalsAlt(0)
    lmm.setNumIntervals0(NumIntervalsDelta0)

    for iA in range(len(Asnps)):
        # add SNP design
        lmm.setSNPcoldesign(Asnps[iA])
        lmm.process()
        pv[iA, :] = lmm.getPv()[0]
    return lmm, pv


def qtl_test_interaction_lmm_kronecker(
    snps,
    phenos,
    covs=None,
    Acovs=None,
    Asnps1=None,
    Asnps0=None,
    K1r=None,
    K1c=None,
    K2r=None,
    K2c=None,
    NumIntervalsDelta0=100,
    NumIntervalsDeltaAlt=100,
    searchDelta=False,
    return_lmm=False,
):
    r"""
    Wrapper function for single-variant multi-trait interaction test

    Args:
        snps (ndarray):
            (`N`, `S`) ndarray of `S` SNPs for `N` individuals.
        pheno (ndarray):
            (`N`, `P`) ndarray of `P` phenotype sfor `N` individuals.
            If phenotypes have missing values, then the subset of
            individuals used for each phenotype column will be subsetted.
        covs (ndarray or list of ndarrays, optional):
            sample design matrices of the different fixed effect terms.
            See `limix.qtl.qtl_test_lmm_kronecker` for more info.
        Acovs (ndarray or list of ndarrays, optional):
            trait design matrices of the different fixed effect terms.
            See `limix.qtl.qtl_test_lmm_kronecker` for more info.
        Asnps1 (ndarray or list of ndarray, optional):
            (`K`, `P`) trait design of rank `K` of the snp fixed effect
            on the `P` traits for the alternative model.
            ``Asnps1`` is tested against the null trait design ``Asnps0``.
            Multiple designs can be tested vs the null trait design
            If a list of ndarray trait designs is specified for ``Asnps1``,
            the multiple designs are tested vs the null trait design
            ``Asnps0`` one-by-one.
            By default ``Asnps1`` is a (`P`, `P`) identity ndarray.
        Asnps0 (ndarray, optional):
            (`K0`, `P`) trait design of rank `K` of the snp fixed effect
            on the `P` traits for the null model.
            By default ``Asnps0`` is (`1`, `P`) ndrray of ones.
        K1r (ndarray, optional):
            (`N`, `N`) ndarray of LMM-covariance/kinship coefficients.
            If not provided, ``K1r = numpy.dot(snps, snps.T)`` is assumed.
        K1c (ndarray, optional):
            (`P`, `P`) ndarray of LMM-covariance/kinship coefficients.
            If either `K1c` or `K2c` is not provided, K1c is estimated by
            (restricted) maximum likelihood  from the model without the snp effect.
        K2r (ndarray, optional):
            (`N`, `N`) ndarray of LMM-covariance/kinship coefficients.
            By default ``K2r`` is a (`N`, `N`) identity ndarray.
        K2c (ndarray, optional):
            (`P`, `P`) ndarray of residual covariances.
            If either `K1c` or `K2c` is not provided, K1c is estimated by
            (restricted) maximum likelihood  from the model without the snp effect.
        NumIntervalsDelta0 (int, optional):
            number of steps for delta optimization on the null model.
            By default ``NumIntervalsDelta0`` is 100.
        NumIntervalsDeltaAlt (int, optional):
            number of steps for delta optimization on the alternative model.
            Requires ``searchDelta=True`` to have an effect.
        searchDelta (bool, optional):
            if True, delta optimization on the alternative model is carried out.
            By default ``searchDelta`` is False.

    Returns:
        (tuple): tuple containing:
            - **pv** (*ndarray*):
              P-values of the interaction test (``Asnps1`` vs ``Asnps0``)
              If ``Asnps1`` is a list of trait designs,
              `pv` contains the P values for the different trait designs (as rows).
            - **pv0** (*ndarray*):
              P-values of the association test ``Asnps0`` vs ``0``.
            - **pvAlt** (*ndarray*):
              P-values of the association test ``Asnps1`` vs ``0``.
              If ``Asnps1`` is a list of trait designs,
              `pvAlt` contains the P values for the different
              trait designs (as rows).

    Examples
    --------

        Example of how to perform a specific, common and any effect test.
        For more information see `limix_tutorials.`_

        .. _limix_tutorials.: http://nbviewer.jupyter.org/github/limix/limix-tutorials/blob/master/index.ipynb

        .. doctest::

            >>> from numpy.random import RandomState
            >>> from numpy import dot, eye, ones
            >>> from limix.qtl import qtl_test_interaction_lmm_kronecker
            >>> random = RandomState(1)
            >>>
            >>> N = 100
            >>> S = 1000
            >>> P = 2
            >>>
            >>> snps = (random.rand(N, S) < 0.2).astype(float)
            >>> pheno = random.randn(N, P)
            >>> W = random.randn(N, 10)
            >>> kinship = dot(W, W.T) / float(10)
            >>> kinship+= 1e-4 * eye(N)
            >>>
            >>> #interaction effect test
            >>> Asnps0 = ones((1,P))
            >>> Asnps1 = eye(P)
            >>> pv, pv0, pvAlt = qtl_test_interaction_lmm_kronecker(snps, pheno,
            ...                                                     K1r=kinship,
            ...                                                     Asnps0=Asnps0,
            ...                                                     Asnps1=Asnps1)
            >>> #interaction P value
            >>> print(pv.shape)
            (1, 1000)
            >>> print(pv[:,:4])
            [[0.9347 0.7744 0.2678 0.2893]]
            >>>
            >>> #common effect P value
            >>> print(pv0.shape)
            (1, 1000)
            >>> print(pv0[:,:4])
            [[0.8722 0.6509 0.3593 0.6816]]
            >>>
            >>> #any effect P value
            >>> print(pvAlt.shape)
            (1, 1000)
            >>> print(pvAlt[:,:4])
            [[0.9838 0.8664 0.3555 0.5244]]
    """
    import limix_legacy.deprecated
    import limix_legacy.deprecated as dlimix_legacy

    S = snps.shape[1]
    # 0. checks
    N = phenos.shape[0]
    P = phenos.shape[1]

    if K1r is None:
        K1r = np.dot(snps, snps.T)
    else:
        assert K1r.shape[0] == N, "K1r: dimensions dismatch"
        assert K1r.shape[1] == N, "K1r: dimensions dismatch"

    if K2r is None:
        K2r = np.eye(N)
    else:
        assert K2r.shape[0] == N, "K2r: dimensions dismatch"
        assert K2r.shape[1] == N, "K2r: dimensions dismatch"

    covs, Acovs = _updateKronCovs(covs, Acovs, N, P)

    # Asnps can be several designs
    if Asnps0 is None:
        Asnps0 = [np.ones([1, P])]
    if Asnps1 is None:
        Asnps1 = [np.eye([P])]
    if not isinstance(Asnps0, list):
        Asnps0 = [Asnps0]
    if not isinstance(Asnps1, list):
        Asnps1 = [Asnps1]
    assert (len(Asnps0) == 1) and (
        len(Asnps1) > 0
    ), "need at least one Snp design matrix for null and alt model"

    # one row per column design matrix
    pv = np.zeros((len(Asnps1), snps.shape[1]))
    lrt = np.zeros((len(Asnps1), snps.shape[1]))
    pvAlt = np.zeros((len(Asnps1), snps.shape[1]))
    lrtAlt = np.zeros((len(Asnps1), snps.shape[1]))

    # 1. run GP model to infer suitable covariance structure
    if K1c is None or K2c is None:
        vc = _estimateKronCovariances(
            phenos=phenos,
            K1r=K1r,
            K2r=K2r,
            K1c=K1c,
            K2c=K2c,
            covs=covs,
            Acovs=Acovs,
            trait_covar_type="freeform",
            verbose=False,
        )
        K1c = vc.getTraitCovar(0)
        K2c = vc.getTraitCovar(1)
    else:
        assert K1c.shape[0] == P, "K1c: dimensions dismatch"
        assert K1c.shape[1] == P, "K1c: dimensions dismatch"
        assert K2c.shape[0] == P, "K2c: dimensions dismatch"
        assert K2c.shape[1] == P, "K2c: dimensions dismatch"

    # 2. run kroneckerLMM for null model
    lmm = limix_legacy.deprecated.CKroneckerLMM()
    lmm.setK1r(K1r)
    lmm.setK1c(K1c)
    lmm.setK2r(K2r)
    lmm.setK2c(K2c)
    lmm.setSNPs(snps)
    # add covariates
    for ic in range(len(Acovs)):
        lmm.addCovariates(covs[ic], Acovs[ic])
    lmm.setPheno(phenos)

    # delta serch on alt. model?
    if searchDelta:
        lmm.setNumIntervalsAlt(NumIntervalsDeltaAlt)
        lmm.setNumIntervals0_inter(NumIntervalsDeltaAlt)
    else:
        lmm.setNumIntervalsAlt(0)
        lmm.setNumIntervals0_inter(0)

    lmm.setNumIntervals0(NumIntervalsDelta0)
    # add SNP design
    lmm.setSNPcoldesign0_inter(Asnps0[0])
    for iA in range(len(Asnps1)):
        lmm.setSNPcoldesign(Asnps1[iA])
        lmm.process()

        pvAlt[iA, :] = lmm.getPv()[0]
        pv[iA, :] = lmm.getPv()[1]
        pv0 = lmm.getPv()[2][np.newaxis, :]
    if return_lmm:
        return pv, pv0, pvAlt, lmm
    else:
        return pv, pv0, pvAlt


def qtl_test_interaction_lmm(
    snps, pheno, Inter, Inter0=None, covs=None, K=None, test="lrt"
):
    r"""
    Wrapper function for single-variant interaction test using LMMs.

    Args:
        snps (ndarray):
            (`N`, `S`) ndarray of `S` SNPs for `N` individuals.
        pheno (ndarray):
            (`N`, `P`) ndarray of `P` phenotype sfor `N` individuals.
            If phenotypes have missing values, then the subset of
            individuals used for each phenotype column will be subsetted.
        Inter (ndarray, optional):
            (`N`, `I`) ndarray of `I` interaction variables to be tested
            for `N` individuals.
        Inter0 (ndarray, optional):
            (`N`, `I0`) ndarray of `I0` interaction variables to be included
            for `N` individuals to be included both in the alternative and
            in the null model.
            By default `Inter0` is a (`N`, `1`) array of ones.
        covs (ndarray, optional):
            (`N`, `D`) ndarray of `D` covariates for `N` individuals.
            By default, ``covs`` is a (`N`, `1`) array of ones.
        K (ndarray, optional):
            (`N`, `N`) ndarray of LMM-covariance/kinship coefficients (optional)
            If not provided, then standard linear regression is considered.
        test ({'lrt', 'f'}, optional):
            test statistic.
            'lrt' for likelihood ratio test (default) or 'f' for F-test.

    Returns:
        :class:`limix_legacy.deprecated.CInteractLMM()`:
            limix CInteractLMM object

    Examples
    --------

        .. doctest::

            >>> from numpy.random import RandomState
            >>> from numpy import dot, eye, ones, concatenate
            >>> from limix.qtl import qtl_test_interaction_lmm
            >>> random = RandomState(1)
            >>>
            >>> N = 100
            >>> S = 1000
            >>>
            >>> snps = (random.rand(N, S) < 0.2).astype(float)
            >>> pheno = random.randn(N, 1)
            >>> W = random.randn(N, 10)
            >>> kinship = dot(W, W.T) / float(10)
            >>> kinship+= 1e-4 * eye(N)
            >>>
            >>> #environments to include in null and alt models
            >>> Inter0 = ones((N,1))
            >>>
            >>> #environments to include in alt model
            >>> Inter = random.randn(N,2)
            >>>
            >>> #environments are also included as fixed effect covariates
            >>> mean = ones([N, 1])
            >>> covs = concatenate([mean, Inter], 1)
            >>>
            >>> #interaction test
            >>> lmi = qtl_test_interaction_lmm(snps, pheno, Inter=Inter,
            ...                                Inter0=Inter0)
            >>> pvi = lmi.getPv()
            >>>
            >>> print(pvi.shape)
            (1, 1000)
            >>> print(pvi[:,:4])
            [[0.8179 0.2185 0.3338 0.077 ]]
    """

    import limix_legacy.deprecated
    import limix_legacy.deprecated as dlimix_legacy

    N = snps.shape[0]
    if covs is None:
        covs = np.ones((N, 1))
    if K is None:
        K = np.eye(N)
    if Inter0 is None:
        Inter0 = np.ones([N, 1])
    assert (
        pheno.shape[0] == N
        and K.shape[0] == N
        and K.shape[1] == N
        and covs.shape[0] == N
        and Inter0.shape[0] == N
        and Inter.shape[0] == N
    ), "shapes missmatch"
    lmi = limix_legacy.deprecated.CInteractLMM()
    lmi.setK(K)
    lmi.setSNPs(snps)
    lmi.setPheno(pheno)
    lmi.setCovs(covs)
    lmi.setInter0(Inter0)
    lmi.setInter(Inter)
    if test == "lrt":
        lmi.setTestStatistics(lmi.TEST_LRT)
    elif test == "f":
        lmi.setTestStatistics(lmi.TEST_F)
    else:
        print(test)
        raise NotImplementedError("only f or lrt are implemented")
    lmi.process()
    return lmi


""" MULTI LOCUS MODEL """


def forward_lmm(
    snps,
    pheno,
    K=None,
    covs=None,
    qvalues=False,
    threshold=5e-8,
    maxiter=2,
    test="lrt",
    verbose=None,
    **kw_args
):
    r"""
    Wrapper function for univariate single-variant test with forward selection

    Args:
        snps (ndarray):
            (`N`, `S`) ndarray of `S` SNPs for `N` individuals.
        pheno (ndarray):
            (`N`, `P`) ndarray of `P` phenotype sfor `N` individuals.
            If phenotypes have missing values, then the subset of
            individuals used for each phenotype column will be subsetted.
        K (ndarray, optional):
            (`N`, `N`) ndarray of LMM-covariance/kinship coefficients (optional)
            If not provided, then standard linear regression is considered.
        covs (ndarray, optional):
            (`N`, `D`) ndarray of `D` covariates for `N` individuals.
            By default, ``covs`` is a (`N`, `1`) array of ones.
        qvalues (bool, optional):
            if True, ``threshold`` is set to Storey qvalues to control FDR.
            (default is False)
        threshold (float, optional):
            P-value thrashold for inclusion in forward selection (default 5e-8).
            If ``qvalues=True``, the threshold is set to to Storey qvalues.
        maxiter (int, optional):
            maximum number of interaction scans. First scan is
            without inclusion, so maxiter-1 inclusions can be performed.
            (default 2)
        test ({'lrt', 'f'}, optional):
            test statistic.
            'lrt' for likelihood ratio test (default) or 'f' for F-test.
        verbose (bool, optional):
            print verbose output. (default is False)

    Returns:
        (tuple): tuple containing:
            - **lmm** (*:class:`limix.qtl.LMM`*): LIMIX LMM object
            - **res** (*dict*): {**iadded**, **pvadded**, **pvall**}.
              **iadded** is an ndarray of indices of SNPs included in order
              of inclusion, **pvadded** is an ndarray of Pvalues obtained by
              the included SNPs in iteration and **pvall** is a  (`Nadded`, `S`)
              ndarray of Pvalues for all iterations

    Examples
    --------

        .. doctest::

            >>> from numpy.random import RandomState
            >>> from numpy import dot, eye, ones
            >>> from limix.qtl import forward_lmm
            >>> random = RandomState(1)
            >>>
            >>> N = 100
            >>> S = 1000
            >>>
            >>> snps = (random.rand(N, S) < 0.2).astype(float)
            >>> pheno = random.randn(N, 1)
            >>> W = random.randn(N, 10)
            >>> kinship = dot(W, W.T) / float(10)
            >>> kinship+= 1e-4 * eye(N)
            >>>
            >>> #forward lmm
            >>> lmm, res = forward_lmm(snps, pheno, K=kinship, threshold=1.)
            >>>
            >>> print(res['pvall'].shape)
            (2, 1000)
            >>> print(res['pvall'][:,:4])
            [[0.8571 0.4668 0.5872 0.5589]
             [0.77   0.4226 0.6165 0.8727]]
    """

    import limix_legacy.deprecated
    import limix_legacy.deprecated as dlimix_legacy

    verbose = dlimix_legacy.getVerbose(verbose)

    if K is None:
        K = np.eye(snps.shape[0])
    if covs is None:
        covs = np.ones((snps.shape[0], 1))
    # assert single trait
    assert pheno.shape[1] == 1, "forward_lmm only supports single phenotypes"

    lm = qtl_test_lmm(snps, pheno, K=K, covs=covs, test=test, **kw_args)
    pvall = []
    pv = lm.getPv().ravel()
    # hack to avoid issues with degenerate pv
    pv[sp.isnan(pv)] = 1
    pvall.append(pv)
    imin = pv.argmin()
    niter = 1
    # start stuff
    iadded = []
    pvadded = []
    qvadded = []
    if qvalues:
        assert (
            pv.shape[0] == 1
        ), "This is untested with the fdr package. pv.shape[0]==1 failed"
        qvall = []
        qv = qvalues(pv)
        qvall.append(qv)
        score = qv.min()
    else:
        score = pv.min()
    while (score < threshold) and niter < maxiter:
        t0 = time.time()
        iadded.append(imin)
        pvadded.append(pv[imin])
        if qvalues:
            qvadded.append(qv[0, imin])
        covs = np.concatenate((covs, snps[:, imin : (imin + 1)]), 1)
        lm.setCovs(covs)
        lm.process()
        pv = lm.getPv().ravel()
        pv[sp.isnan(pv)] = 1
        pvall.append(pv)
        imin = pv.argmin()
        if qvalues:
            qv = qvalues(pv)
            qvall[niter : niter + 1, :] = qv
            score = qv.min()
        else:
            score = pv.min()
        t1 = time.time()
        if verbose:
            print(("finished GWAS testing in %.2f seconds" % (t1 - t0)))
        niter = niter + 1
    RV = {}
    RV["iadded"] = iadded
    RV["pvadded"] = pvadded
    RV["pvall"] = np.array(pvall)
    if qvalues:
        RV["qvall"] = np.array(qvall)
        RV["qvadded"] = qvadded
    return lm, RV


# TOOD: use **kw_args to forward params.. see below
def forward_lmm_kronecker(
    snps,
    phenos,
    Asnps=None,
    Acond=None,
    K1r=None,
    K1c=None,
    K2r=None,
    K2c=None,
    covs=None,
    Acovs=None,
    threshold=5e-8,
    maxiter=2,
    qvalues=False,
    update_covariances=False,
    verbose=None,
    **kw_args
):
    """
    Kronecker fixed effects test with forward selection

    Args:
        snps:   [N x S] np.array of S SNPs for N individuals (test SNPs)
        pheno:  [N x P] np.array of 1 phenotype for N individuals
        K:      [N x N] np.array of LMM-covariance/kinship koefficients (optional)
                        If not provided, then linear regression analysis is performed
        covs:   [N x D] np.array of D covariates for N individuals
        threshold:      (float) P-value thrashold for inclusion in forward selection (default 5e-8)
        maxiter:        (int) maximum number of interaction scans. First scan is
                        without inclusion, so maxiter-1 inclusions can be performed. (default 2)
        qvalues:        Use q-value threshold and return q-values in addition (default False)
        update_covar:   Boolean indicator if covariances should be re-estimated after each forward step (default False)

    Returns:
        lm:             lmix LMMi object
        resultStruct with elements:
            iadded:         array of indices of SNPs included in order of inclusion
            pvadded:        array of Pvalues obtained by the included SNPs in iteration
                            before inclusion
            pvall:     [Nadded x S] np.array of Pvalues for all iterations.
        Optional:      corresponding q-values
            qvadded
            qvall
    """
    import limix_legacy.deprecated
    import limix_legacy.deprecated as dlimix_legacy

    verbose = dlimix_legacy.getVerbose(verbose)
    # 0. checks
    N = phenos.shape[0]
    P = phenos.shape[1]

    if K1r is None:
        K1r = np.dot(snps, snps.T)
    else:
        assert K1r.shape[0] == N, "K1r: dimensions dismatch"
        assert K1r.shape[1] == N, "K1r: dimensions dismatch"

    if K2r is None:
        K2r = np.eye(N)
    else:
        assert K2r.shape[0] == N, "K2r: dimensions dismatch"
        assert K2r.shape[1] == N, "K2r: dimensions dismatch"

    covs, Acovs = _updateKronCovs(covs, Acovs, N, P)

    if Asnps is None:
        Asnps = [np.ones([1, P])]
    if not isinstance(Asnps, list):
        Asnps = [Asnps]
    assert len(Asnps) > 0, "need at least one Snp design matrix"

    if Acond is None:
        Acond = Asnps
    if not isinstance(Acond, list):
        Acond = [Acond]
    assert len(Acond) > 0, "need at least one Snp design matrix"

    # 1. run GP model to infer suitable covariance structure
    if K1c is None or K2c is None:
        vc = _estimateKronCovariances(
            phenos=phenos,
            K1r=K1r,
            K2r=K2r,
            K1c=K1c,
            K2c=K2c,
            covs=covs,
            Acovs=Acovs,
            **kw_args
        )
        K1c = vc.getTraitCovar(0)
        K2c = vc.getTraitCovar(1)
    else:
        vc = None
        assert K1c.shape[0] == P, "K1c: dimensions dismatch"
        assert K1c.shape[1] == P, "K1c: dimensions dismatch"
        assert K2c.shape[0] == P, "K2c: dimensions dismatch"
        assert K2c.shape[1] == P, "K2c: dimensions dismatch"
    t0 = time.time()
    lm, pv = qtl_test_lmm_kronecker(
        snps=snps,
        phenos=phenos,
        Asnps=Asnps,
        K1r=K1r,
        K2r=K2r,
        K1c=K1c,
        K2c=K2c,
        covs=covs,
        Acovs=Acovs,
    )

    # get pv
    # start stuff
    iadded = []
    pvadded = []
    qvadded = []
    time_el = []
    pvall = []
    qvall = None
    t1 = time.time()
    if verbose:
        print(("finished GWAS testing in %.2f seconds" % (t1 - t0)))
    time_el.append(t1 - t0)
    pvall.append(pv)
    imin = np.unravel_index(pv.argmin(), pv.shape)
    score = pv[imin].min()
    niter = 1
    if qvalues:
        assert (
            pv.shape[0] == 1
        ), "This is untested with the fdr package. pv.shape[0]==1 failed"
        qvall = []
        qv = qvalues(pv)
        qvall.append(qv)
        score = qv[imin]
    # loop:
    while (score < threshold) and niter < maxiter:
        t0 = time.time()
        pvadded.append(pv[imin])
        iadded.append(imin)
        if qvalues:
            qvadded.append(qv[imin])
        if update_covariances and vc is not None:
            vc.addFixedTerm(snps[:, imin[1] : (imin[1] + 1)], Acond[imin[0]])
            vc.setScales()  # CL: don't know what this does, but findLocalOptima crashes becahuse vc.noisPos=None
            vc.findLocalOptima(fast=True)
            K1c = vc.getTraitCovar(0)
            K2c = vc.getTraitCovar(1)
            lm.setK1c(K1c)
            lm.setK2c(K2c)
        lm.addCovariates(snps[:, imin[1] : (imin[1] + 1)], Acond[imin[0]])
        for i in range(len(Asnps)):
            # add SNP design
            lm.setSNPcoldesign(Asnps[i])
            lm.process()
            pv[i, :] = lm.getPv()[0]
        pvall.append(pv.ravel())
        imin = np.unravel_index(pv.argmin(), pv.shape)
        if qvalues:
            qv = qvalues(pv)
            qvall[niter : niter + 1, :] = qv
            score = qv[imin].min()
        else:
            score = pv[imin].min()
        t1 = time.time()
        if verbose:
            print(("finished GWAS testing in %.2f seconds" % (t1 - t0)))
        time_el.append(t1 - t0)
        niter = niter + 1
    RV = {}
    RV["iadded"] = iadded
    RV["pvadded"] = pvadded
    RV["pvall"] = np.array(pvall)
    RV["time_el"] = time_el
    if qvalues:
        RV["qvall"] = qvall
        RV["qvadded"] = qvadded
    return lm, RV


""" INTERNAL """


def _estimateKronCovariances(
    phenos,
    K1r=None,
    K1c=None,
    K2r=None,
    K2c=None,
    covs=None,
    Acovs=None,
    trait_covar_type="freeform",
    rank=1,
    lambd=None,
    verbose=True,
    init_method="random",
    old_opt=True,
):
    """
    estimates the background covariance model before testing

    Args:
        phenos: [N x P] np.array of P phenotypes for N individuals
        K1r:    [N x N] np.array of LMM-covariance/kinship koefficients (optional)
                        If not provided, then linear regression analysis is performed
        K1c:    [P x P] np.array of LMM-covariance/kinship koefficients (optional)
                        If not provided, then linear regression analysis is performed
        K2r:    [N x N] np.array of LMM-covariance/kinship koefficients (optional)
                        If not provided, then linear regression analysis is performed
        K2c:    [P x P] np.array of LMM-covariance/kinship koefficients (optional)
                        If not provided, then linear regression analysis is performed
        covs:           list of np.arrays holding covariates. Each covs[i] has one corresponding Acovs[i]
        Acovs:          list of np.arrays holding the phenotype design matrices for covariates.
                        Each covs[i] has one corresponding Acovs[i].
        trait_covar_type:     type of covaraince to use. Default 'freeform'. possible values are
                        'freeform': free form optimization,
                        'fixed': use a fixed matrix specified in covar_K0,
                        'diag': optimize a diagonal matrix,
                        'lowrank': optimize a low rank matrix. The rank of the lowrank part is specified in the variable rank,
                        'lowrank_id': optimize a low rank matrix plus the weight of a constant diagonal matrix. The rank of the lowrank part is specified in the variable rank,
                        'lowrank_diag': optimize a low rank matrix plus a free diagonal matrix. The rank of the lowrank part is specified in the variable rank,
                        'block': optimize the weight of a constant P x P block matrix of ones,
                        'block_id': optimize the weight of a constant P x P block matrix of ones plus the weight of a constant diagonal matrix,
                        'block_diag': optimize the weight of a constant P x P block matrix of ones plus a free diagonal matrix,
        rank:           rank of a possible lowrank component (default 1)

    Returns:
        VarianceDecomposition object
    """
    import limix_legacy.modules.varianceDecomposition as VAR

    # from . import varianceDecomposition as VAR
    if verbose:
        print(".. Estimating trait covariances")
    vc = VAR.VarianceDecomposition(phenos)
    if K1r is not None:
        vc.addRandomEffect(K1r, trait_covar_type="freeform", rank=rank)
    if K2r is not None:
        # TODO: fix this; forces second term to be the noise covariance
        vc.addRandomEffect(is_noise=True, K=K2r, trait_covar_type="freeform", rank=rank)
    for ic in range(len(Acovs)):
        vc.addFixedEffect(covs[ic], Acovs[ic])
    start = time.time()
    if old_opt:
        conv = vc.optimize(fast=True)
    elif lambd is not None:
        conv = vc.optimize(init_method=init_method, verbose=verbose, lambd=lambd)
    else:
        conv = vc.optimize(init_method=init_method, verbose=verbose)
    assert conv, "Variance Decomposition has not converged"
    time_el = time.time() - start
    if verbose:
        print(("Done in %.2f s" % time_el))
    return vc


def _updateKronCovs(covs, Acovs, N, P):
    """
    make sure that covs and Acovs are lists
    """
    if (covs is None) and (Acovs is None):
        covs = [np.ones([N, 1])]
        Acovs = [np.eye(P)]

    if Acovs is None or covs is None:
        raise Exception("Either Acovs or covs is None, while the other isn't")

    if (not isinstance(Acovs, list)) and (not isinstance(covs, list)):
        Acovs = [Acovs]
        covs = [covs]
    if (
        (not isinstance(covs, list))
        or (not isinstance(Acovs, list))
        or (len(covs) != len(Acovs))
    ):
        raise Exception(
            "Either Acovs or covs is not a list or they missmatch in length"
        )
    return covs, Acovs


# TODO: we need to fix. THis does not work as interact_GxE is not existing
# I vote we also use **kw_args to forward parameters to interact_Gxe?
def qtl_test_interaction_GxG(pheno, snps1, snps2=None, K=None, covs=None, test="lrt"):
    """
    Epistasis test between two sets of SNPs

    Args:
        pheno:  [N x 1] np.array of 1 phenotype for N individuals
        snps1:  [N x S1] np.array of S1 SNPs for N individuals
        snps2:  [N x S2] np.array of S2 SNPs for N individuals
        K:      [N x N] np.array of LMM-covariance/kinship koefficients (optional)
                        If not provided, then linear regression analysis is performed
        covs:   [N x D] np.array of D covariates for N individuals
        test:    'lrt' for likelihood ratio test (default) or 'f' for F-test

    Returns:
        pv:     [S2 x S1] np.array of P values for epistasis tests beten all SNPs in
                snps1 and snps2
    """
    if K is None:
        K = np.eye(N)
    N = snps1.shape[0]
    if snps2 is None:
        snps2 = snps1
    return qtl_test_interaction_GxE_1dof(
        snps=snps1, pheno=pheno, env=snps2, covs=covs, K=K, test=test
    )


def qtl_test_interaction_GxE_1dof(
    snps, pheno, env, K=None, covs=None, test="lrt", verbose=None
):
    """
    Univariate GxE fixed effects interaction linear mixed model test for all
    pairs of SNPs and environmental variables.

    Args:
        snps:   [N x S] np.array of S SNPs for N individuals
        pheno:  [N x 1] np.array of 1 phenotype for N individuals
        env:    [N x E] np.array of E environmental variables for N individuals
        K:      [N x N] np.array of LMM-covariance/kinship koefficients (optional)
                        If not provided, then linear regression analysis is performed
        covs:   [N x D] np.array of D covariates for N individuals
        test:    'lrt' for likelihood ratio test (default) or 'f' for F-test
        verbose: print verbose output? (False)

    Returns:
        pv:     [E x S] np.array of P values for interaction tests between all
                E environmental variables and all S SNPs
    """
    import limix_legacy.deprecated as dlimix_legacy

    verbose = dlimix_legacy.getVerbose(verbose)
    N = snps.shape[0]
    if K is None:
        K = np.eye(N)
    if covs is None:
        covs = np.ones((N, 1))
    assert (
        env.shape[0] == N
        and pheno.shape[0] == N
        and K.shape[0] == N
        and K.shape[1] == N
        and covs.shape[0] == N
    ), "shapes missmatch"
    Inter0 = np.ones((N, 1))
    pv = np.zeros((env.shape[1], snps.shape[1]))
    if verbose:
        print(
            (
                "starting %i interaction scans for %i SNPs each."
                % (env.shape[1], snps.shape[1])
            )
        )
    t0 = time.time()
    for i in range(env.shape[1]):
        t0_i = time.time()
        cov_i = np.concatenate((covs, env[:, i : (i + 1)]), 1)
        lm_i = qtl_test_interaction_lmm(
            snps=snps,
            pheno=pheno,
            covs=cov_i,
            Inter=env[:, i : (i + 1)],
            Inter0=Inter0,
            test=test,
        )
        pv[i, :] = lm_i.getPv()[0, :]
        t1_i = time.time()
        if verbose:
            print(
                (
                    "Finished %i out of %i interaction scans in %.2f seconds."
                    % ((i + 1), env.shape[1], (t1_i - t0_i))
                )
            )
    t1 = time.time()
    print(
        (
            "-----------------------------------------------------------\nFinished all %i interaction scans in %.2f seconds."
            % (env.shape[1], (t1 - t0))
        )
    )
    return pv


def phenSpecificEffects(snps, pheno1, pheno2, K=None, covs=None, test="lrt"):
    """
    Univariate fixed effects interaction test for phenotype specific SNP effects

    Args:
        snps:   [N x S] np.array of S SNPs for N individuals (test SNPs)
        pheno1: [N x 1] np.array of 1 phenotype for N individuals
        pheno2: [N x 1] np.array of 1 phenotype for N individuals
        K:      [N x N] np.array of LMM-covariance/kinship koefficients (optional)
                        If not provided, then linear regression analysis is performed
        covs:   [N x D] np.array of D covariates for N individuals
        test:    'lrt' for likelihood ratio test (default) or 'f' for F-test

    Returns:
        limix LMM object
    """
    N = snps.shape[0]
    if K is None:
        K = np.eye(N)
    assert (
        pheno1.shape[1] == pheno2.shape[1]
    ), "Only consider equal number of phenotype dimensions"
    if covs is None:
        covs = np.ones(N, 1)
    assert (
        pheno1.shape[1] == 1
        and pheno2.shape[1] == 1
        and pheno1.shape[0] == N
        and pheno2.shape[0] == N
        and K.shape[0] == N
        and K.shape[1] == N
        and covs.shape[0] == N
    ), "shapes missmatch"
    Inter = np.zeros((N * 2, 1))
    Inter[0:N, 0] = 1
    Inter0 = np.ones((N * 2, 1))
    Yinter = np.concatenate((pheno1, pheno2), 0)
    Xinter = np.tile(snps, (2, 1))
    Covitner = np.tile(covs(2, 1))
    lm = qtl_test_interaction_lmm(
        snps=Xinter, pheno=Yinter, covs=Covinter, Inter=Inter, Inter0=Inter0, test=test
    )
    return lm
