import numpy as np
import scipy.stats as st
import scipy as sp
import scipy.linalg as la
import pandas as pd
import time
from struct_lmm import StructLMM
from limix_core.gp import GP2KronSum
from limix_core.gp import GP2KronSumLR
from limix_core.covar import FreeFormCov


def add_jitter(S_R):
    assert S_R.min() > -1e-6, "LMM-covariance is not sdp!"
    RV = S_R + sp.maximum(1e-4 - S_R.min(), 0)
    return RV


class GWAS_LMM:
    """
    Wrapper function for univariate single-variant association testing
    using variants of the linear mixed model.

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

    def __init__(
        self,
        pheno,
        covs=None,
        R=None,
        eigh_R=None,
        W_R=None,
        inter=None,
        inter0=None,
        verbose=False,
    ):
        from limix_lmm.lmm import LMM
        from limix_lmm.lmm_core import LMMCore

        self.verbose = None

        if covs is None:
            covs = sp.ones([pheno.shape[0], 1])

        # case 1: linear model
        if W_R is None and eigh_R is None and R is None:
            if verbose:
                print("Model: lm")
            self.gp = None
            Kiy_fun = None

        # case 2: low-rank linear model
        elif W_R is not None:
            if verbose:
                print("Model: low-rank lmm")
            self.gp = GP2KronSumLR(
                Y=pheno, Cn=FreeFormCov(1), G=W_R, F=covs, A=sp.ones((1, 1))
            )
            self.gp.covar.Cr.setCovariance(sp.var(pheno) * sp.ones((1, 1)))
            self.gp.covar.Cn.setCovariance(sp.var(pheno) * sp.ones((1, 1)))
            info_opt = self.gp.optimize(verbose=verbose)
            Kiy_fun = self.gp.covar.solve

        # case 3: full-rank linear model
        else:
            if verbose:
                print("Model: lmm")
            if eigh_R is None:
                eigh_R = la.eigh(R)
            S_R, U_R = eigh_R
            add_jitter(S_R)
            self.gp = GP2KronSum(
                Y=pheno,
                Cg=FreeFormCov(1),
                Cn=FreeFormCov(1),
                S_R=S_R,
                U_R=U_R,
                F=covs,
                A=sp.ones((1, 1)),
            )
            self.gp.covar.Cr.setCovariance(0.5 * sp.var(pheno) * sp.ones((1, 1)))
            self.gp.covar.Cn.setCovariance(0.5 * sp.var(pheno) * sp.ones((1, 1)))
            info_opt = self.gp.optimize(verbose=verbose)
            Kiy_fun = self.gp.covar.solve

        if inter is None:
            self.lmm = LMM(pheno, covs, Kiy_fun)
            self.inter1 = None
            self.inter0 = None
        else:
            self.lmm = LMMCore(pheno, covs, Kiy_fun)
            if inter0 is None:
                inter0 = sp.ones([pheno.shape[0], 1])
            if (inter0 == 1).sum():
                self.lmm0 = LMM(pheno, covs, Kiy_fun)
            else:
                self.lmm0 = LMMCore(pheno, covs, Kiy_fun)
            self.inter1 = sp.concatenate([inter0, inter], 1)
            self.inter0 = inter0

    def process(self, snps, return_ste=False, return_lrt=False):
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
        if self.inter1 is None:

            self.lmm.process(snps)
            RV = {}
            RV["pv"] = self.lmm.getPv()
            RV["beta"] = self.lmm.getBetaSNP()
            if return_ste:
                RV["beta_ste"] = self.lmm.getBetaSNPste()
            if return_lrt:
                RV["lrt"] = self.lmm.getLRT()

        else:

            self.lmm.process(snps, self.inter1)
            if (self.inter0 == 1).sum():
                self.lmm0.process(snps)
            else:
                self.lmm0.process(snps, self.inter0)

            # compute pv
            lrt1 = self.lmm.getLRT()
            lrt0 = self.lmm0.getLRT()
            lrt = lrt1 - lrt0
            pv = st.chi2(self.inter1.shape[1] - self.inter0.shape[1]).sf(lrt)

            RV = {}
            RV["pv1"] = self.lmm.getPv()
            RV["pv0"] = self.lmm0.getPv()
            RV["pv"] = pv
            if (self.inter0 == 1).sum():
                RV["beta0"] = self.lmm0.getBetaSNP()
                if return_ste:
                    RV["beta0_ste"] = self.lmm0.getBetaSNPste()
            if return_lrt:
                RV["lrt1"] = lrt1
                RV["lrt0"] = lrt0
                RV["lrt"] = lrt

        return pd.DataFrame(RV)


class GWAS_StructLMM:
    """
    Wrapper function for univariate single-variant association testing
    using variants of the linear mixed model.

    Parameters
    ----------
    pheno : (`N`, 1) ndarray
        phenotype data
    environments : (`N`, `E`) ndarray
        environments data.
    covs : (`N`, `D`) ndarray
        covariate design matrix.
        By default, ``covs`` is a (`N`, `1`) array of ones.
    W_R : (`N`, `R`) ndarray
        ``W_R`` define the lmm-covariance as
        ``R`` = dot(``W_R``, transpose(``W_R``)).
        By detault, ``W_R`` is set to ``environments``.
    tests : list
        Which tests are performed.
        Element list values are ``'inter'`` and ``'assoc'``.
        By default, only the interaction test is considered.
    rhos : list
        for the association test, a list of ``rho`` values must be specified.
        The choice of ``rho`` affects the statistical power of the test
        (for more information see the StructLMM paper).
        By default, ``rho=[0, 0.1**2, 0.2**2, 0.3**2, 0.4**2, 0.5**2, 0.5, 1.]``
    verbose : (bool, optional):
        if True, details such as runtime as displayed.
    """

    def __init__(
        self,
        pheno,
        environments,
        covs=None,
        W_R=None,
        tests=None,
        rhos=None,
        verbose=False,
    ):
        self.verbose = verbose
        self.tests = tests

        if covs is None:
            covs = sp.ones((env.shape[0], 1))
        self.covs = covs

        if rhos is None:
            rhos = [0.0, 0.1 ** 2, 0.2 ** 2, 0.3 ** 2, 0.4 ** 2, 0.5 ** 2, 0.5, 1.0]

        if W_R is None:
            W_R = environments

        if tests is None:
            tests = ["inter"]

        if "inter" in tests:
            self.slmi = StructLMM(pheno, environments, W=environments, rho_list=[0])

        if "assoc" in tests:
            self.slmm = StructLMM(pheno, environments, W=environments, rho_list=rhos)
            self.slmm.fit_null(F=covs, verbose=False)

    def process(self, snps):
        """
        Parameters
        ----------
        snps : (`N`, `S`) ndarray
            genotype data

        Return
        ------
        res : pandas DataFrame
            Results as pandas dataframs
        """
        _pvi = sp.zeros(snps.shape[1])
        _pva = sp.zeros(snps.shape[1])
        for snp in range(snps.shape[1]):
            x = snps[:, [snp]]

            if "inter" in self.tests:
                # interaction test
                covs1 = sp.hstack((self.covs, x))
                self.slmi.fit_null(F=covs1, verbose=False)
                _pvi[snp] = self.slmi.score_2_dof(x)

            if "assoc" in self.tests:
                # association test
                _pva[snp] = self.slmm.score_2_dof(x)

        return pd.DataFrame({"pvi": _pvi, "pva": _pva})


class GWAS_MTLMM:
    """
    Wrapper function for multi-trait single-variant association testing
    using variants of the multi-trait linear mixed model.

    Parameters
    ----------
    pheno : (`N`, `P`) ndarray
        phenotype data
    Asnps : (`P`, `K`) ndarray
         trait design of snp covariance.
         By default, ``Asnps`` is eye(`P`).
    R : (`N`, `N`) ndarray
        LMM-covariance/genetic relatedness matrix.
        If not provided, then standard linear regression is considered.
        Alternatively, its eighenvalue decomposition can be
        provided through ``eigh_R``.
        if ``eigh_R`` is set, this parameter is ignored.
    eigh_R : tuple
        Tuple with `N` ndarray of eigenvalues of `R` and
        (`N`, `N`) ndarray of eigenvectors of ``R``.
    covs : (`N`, `D`) ndarray
        covariate design matrix.
        By default, ``covs`` is a (`N`, `1`) array of ones.
    Acovs : (`P`, `L`) ndarray
        trait design matrices of the different fixed effect terms.
        By default, ``Acovs`` is eye(`P`).
    Asnps0 : (`P`, `K`) ndarray
         trait design of snp covariance in the null model.
         By default, Asnps0 is not considered (i.e., no SNP effect in the null model).
         If specified, then three tests are considered:
         (i) Asnps vs , (ii) Asnps0!=0, (iii) Asnps!=Asnps0
    verbose : (bool, optional):
        if True, details such as runtime as displayed.
    """

    def __init__(
        self,
        pheno,
        Asnps=None,
        R=None,
        eigh_R=None,
        covs=None,
        Acovs=None,
        verbose=None,
        Asnps0=None,
    ):
        self.verbose = None
        from limix_lmm.mtlmm import MTLMM

        if covs is None:
            covs = sp.ones([pheno.shape[0], 1])

        if Acovs is None:
            Acovs = sp.eye(pheno.shape[1])

        # case 1: multi-trait linear model
        assert not (
            eigh_R is None and R is None
        ), "multi-trait linear model not supported"

        # case 2: full-rank multi-trait linear model
        if eigh_R is None:
            eigh_R = la.eigh(R)
        S_R, U_R = eigh_R
        S_R = add_jitter(S_R)
        self.gp = GP2KronSum(
            Y=pheno,
            Cg=FreeFormCov(pheno.shape[1]),
            Cn=FreeFormCov(pheno.shape[1]),
            S_R=eigh_R[0],
            U_R=eigh_R[1],
            F=covs,
            A=Acovs,
        )
        self.gp.covar.Cr.setCovariance(0.5 * sp.cov(pheno.T))
        self.gp.covar.Cn.setCovariance(0.5 * sp.cov(pheno.T))
        info_opt = self.gp.optimize(verbose=verbose)

        self.lmm = MTLMM(pheno, F=covs, A=Acovs, Asnp=Asnps, covar=self.gp.covar)
        if Asnps0 is not None:
            self.lmm0 = MTLMM(pheno, F=covs, A=Acovs, Asnp=Asnps0, covar=self.gp.covar)

        self.Asnps = Asnps
        self.Asnps0 = Asnps0

    def process(self, snps, return_lrt=False):
        """
        Parameters
        ----------
        snps : (`N`, `S`) ndarray
            genotype data
        return_lrt : bool
            if True, return log likelihood ratio tests (default is False)

        Return
        ------
        res : pandas DataFrame
            Results as pandas dataframs
        """
        if self.Asnps0 is None:

            self.lmm.process(snps)
            RV = {}
            RV["pv"] = self.lmm.getPv()
            if return_lrt:
                RV["lrt"] = self.lmm.getLRT()

        else:

            self.lmm.process(snps)
            self.lmm0.process(snps)

            # compute pv
            lrt1 = self.lmm.getLRT()
            lrt0 = self.lmm0.getLRT()
            lrt = lrt1 - lrt0
            pv = st.chi2(self.Asnps.shape[1] - self.Asnps0.shape[1]).sf(lrt)

            RV = {}
            RV["pv1"] = self.lmm.getPv()
            RV["pv0"] = self.lmm0.getPv()
            RV["pv"] = pv
            if return_lrt:
                RV["lrt1"] = lrt1
                RV["lrt0"] = lrt0
                RV["lrt"] = lrt

        return pd.DataFrame(RV)
