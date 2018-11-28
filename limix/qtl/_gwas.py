

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
        from numpy import ones, eye, cov
        from scipy.linalg import eigh
        from limix_core.gp import GP2KronSum
        from limix_core.covar import FreeFormCov

        self.verbose = None
        from limix_lmm.mtlmm import MTLMM

        if covs is None:
            covs = ones([pheno.shape[0], 1])

        if Acovs is None:
            Acovs = eye(pheno.shape[1])

        # case 1: multi-trait linear model
        assert not (
            eigh_R is None and R is None
        ), "multi-trait linear model not supported"

        # case 2: full-rank multi-trait linear model
        if eigh_R is None:
            eigh_R = eigh(R)
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
        self.gp.covar.Cr.setCovariance(0.5 * cov(pheno.T))
        self.gp.covar.Cn.setCovariance(0.5 * cov(pheno.T))
        self.gp.optimize(verbose=verbose)

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
        from pandas import DataFrame
        from scipy.stats import chi2

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
            pv = chi2(self.Asnps.shape[1] - self.Asnps0.shape[1]).sf(lrt)

            RV = {}
            RV["pv1"] = self.lmm.getPv()
            RV["pv0"] = self.lmm0.getPv()
            RV["pv"] = pv
            if return_lrt:
                RV["lrt1"] = lrt1
                RV["lrt0"] = lrt0
                RV["lrt"] = lrt

        return DataFrame(RV)


def add_jitter(S_R):
    from numpy import maximum

    assert S_R.min() > -1e-6, "LMM-covariance is not sdp!"
    RV = S_R + maximum(1e-4 - S_R.min(), 0)
    return RV
