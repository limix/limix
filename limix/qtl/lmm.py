# Copyright(c) 2014, The LIMIX developers (Christoph Lippert, Paolo Francesco Casale, Oliver Stegle)
#
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

import numpy as np
import scipy as sp
import scipy.stats as st
import time

class lmm:
    r"""Basic class for univariate single-variant association testing with LMMs.

    Args:
        snps (ndarray):
            [N, S] ndarray of S SNPs for N individuals.
        pheno (ndarray):
            [N, P] ndarray of P phenotype sfor N individuals.
            If phenotypes have missing values, then the subset of
            individuals used for each phenotype column will be subsetted.
        K (ndarray, optional):
            [N, N] ndarray of LMM-covariance/kinship coefficients (optional)
            If not provided, then standard linear regression is considered.
        covs (ndarray, optional):
            [N, D] ndarray of D covariates for N individuals.
            By default, covs is a [N, 1] vector of ones.
        test ({'lrt', 'f'}, optional):
            test statistic.
            'lrt' for likelihood ratio test (default) or 'f' for F-test.
        NumIntervalsDelta0 (int, optional):
            number of steps for delta optimization on the null model.
            By default `NumIntervalsDelta0` is 100.
        NumIntervalsDeltaAlt (int, optional):
            number of steps for delta optimization on the alternative model.
            Requires `searchDelta`=True to have an effect.
        searchDelta (bool, optional):
            if True, delta optimization on the alternative model is carried out.
            By default `searchDelta` is False.
        verbose (bool, optional):
            if True, details such as runtime as displayed. 

    Example
    -------

        Example of single-variant association testing using a linear mixed model.

        .. doctest::

            >>> from numpy.random import RandomState
            >>> import limix.qtl.lmm
            >>> from numpy import dot
            >>> random = RandomState(1)
            >>>
            >>> N = 100
            >>> S = 1000
            >>>
            >>> # generate data 
            >>> snps = (random.rand(N, S) < 0.2).astype(float)
            >>> pheno = random.randn(N, 1)
            >>> W = random.randn(N, 10)
            >>> kinship = dot(W, W.T) / float(10)
            >>>
            >>> # run single-variant associaiton testing with LMM
            >>> lmm = limix.qtl.lmm(snps, y, kinship)
            >>> pv = lmm.getPv()
            >>> beta = lmm.getBetaSNP()
            >>> beta_ste = lmm.getBetaSNPste()
            >>>
            >>> print lmm.getPv()
            [[ 0.85007143  0.95042782  0.49687393  0.01852559  0.61253936  0.13709002
               0.23412486  0.78590904  0.67982972  0.97131758]]
            >>> print lmm.getBetaSNP()
            [[-0.02113729 -0.00629632 -0.06198703 -0.246229   -0.05527605  0.14154055
              -0.11656051  0.0287622   0.04245761 -0.00370021]]
            >>> print lmm.getBetaSNPste()
            [[ 0.11182137  0.10127666  0.09123567  0.10455841  0.1091437   0.09520354
               0.09796644  0.1058887   0.10287878  0.10290999]]

        As shown in the exmaple below,
        the method can be applied directly on multiple phenotypes.
        Each phenotype is indipendently tested against all variants one-by-one.

        .. doctest::

            >>> random = RandomState(1)
            >>>
            >>> P = 3
            >>> phenos = random.randn(N, P)
            >>>
            >>> lmm = limix.qtl.lmm(snps, phenos, kinship)
            >>> pv = lmm.getPv()
            >>> beta = lmm.getBetaSNP()
            >>> beta_ste = lmm.getBetaSNPste()
            >>>
            >>> print lmm.getPv()
            [[ 0.10443428  0.2211579   0.77353526  0.71451599  0.23743919  0.64876719
               0.62777478  0.17119043  0.81580569  0.12362014]
             [ 0.48500889  0.99372962  0.39939168  0.35453793  0.13346182  0.12371834
               0.73782101  0.40214461  0.23379901  0.28227524]
             [ 0.49768012  0.14925783  0.24932633  0.34892866  0.11523963  0.66897979
               0.62862356  0.62567987  0.99164983  0.36997474]]
    """

    def __init__(self, snps, pheno, K=None, covs=None, test='lrt', NumIntervalsDelta0=100, NumIntervalsDeltaAlt=100, searchDelta=False, verbose=None):
        #create column of 1 for fixed if nothing provide
        try:
            import limix_legacy.deprecated
            import limix_legacy.deprecated as dlimix_legacy
        except ImportError:
            print("Please, install limix-legacy to use this functionality.")

        if len(pheno.shape)==1:
            pheno = pheno[:,sp.newaxis]

        self.verbose = dlimix_legacy.getVerbose(verbose)
        self.snps = snps
        self.pheno = pheno
        self.K = K
        self.covs = covs
        self.test = test
        self.NumIntervalsDelta0 = NumIntervalsDelta0
        self.NumIntervalsDeltaAlt = NumIntervalsDeltaAlt
        self.searchDelta = searchDelta
        self.verbose = verbose
        self.N       = self.pheno.shape[0]
        self.P       = self.pheno.shape[1]
        self.Iok     = ~(np.isnan(self.pheno).any(axis=1))
        if self.K is None:
            self.searchDelta=False
            self.K = np.eye(self.snps.shape[0])
        if self.covs is None:
            self.covs = np.ones((self.snps.shape[0],1))

        self._lmm = None
        #run
        self.verbose = verbose
        self.process()

    def process(self):

        try:
            import limix_legacy.deprecated
            import limix_legacy.deprecated as dlimix_legacy
        except ImportError:
            print("Please, install limix-legacy to use this functionality.")

        t0 = time.time()
        if self._lmm is None:
            self._lmm = limix_legacy.deprecated.CLMM()
            self._lmm.setK(self.K)
            self._lmm.setSNPs(self.snps)
            self._lmm.setPheno(self.pheno)
            self._lmm.setCovs(self.covs)
            if self.test=='lrt':
                self._lmm.setTestStatistics(self._lmm.TEST_LRT)
            elif self.test=='f':
                self._lmm.setTestStatistics(self._lmm.TEST_F)
            else:
                print((self.test))
                raise NotImplementedError("only f and lrt are implemented")
            #set number of delta grid optimizations?
            self._lmm.setNumIntervals0(self.NumIntervalsDelta0)
            if self.searchDelta:
                self._lmm.setNumIntervalsAlt(self.NumIntervalsDeltaAlt)
            else:
                self._lmm.setNumIntervalsAlt(0)

        if not np.isnan(self.pheno).any():
            #process
            self._lmm.process()
            self.pvalues = self._lmm.getPv()
            self.beta_snp = self._lmm.getBetaSNP()
            self.beta_ste = self._lmm.getBetaSNPste()
            self.ldelta_0 = self._lmm.getLdelta0()
            self.ldelta_alt = self._lmm.getLdeltaAlt()
            self.NLL_0 = self._lmm.getNLL0()
            self.NLL_alt = self._lmm.getNLLAlt()
        else:
            if self._lmm is not None:
                raise Exception('cannot reuse a CLMM object if missing variables are present')
            else:
                self._lmm = limix_legacy.deprecated.CLMM()
            #test all phenotypes separately
            self.pvalues = np.zeros((self.phenos.shape[1],self.snps.shape[1]))
            self.beta_snp = np.zeros((self.phenos.shape[1],self.snps.shape[1]))
            self.beta_ste = np.zeros((self.phenos.shape[1],self.snps.shape[1]))
            self.ldelta_0 = np.zeros((self.phenos.shape[1],self.snps.shape[1]))
            self.ldelta_alt = np.zeros((self.phenos.shape[1],self.snps.shape[1]))
            self.NLL_0 = np.zeros((self.phenos.shape[1],self.snps.shape[1]))
            self.NLL_alt = np.zeros((self.phenos.shape[1],self.snps.shape[1]))
            self.test_statistics = np.zeros((self.phenos.shape[1],self.snps.shape[1]))
            for ip in np.arange(self.phenos.shape[1]):
                pheno_ = self.phenos[:,ip]
                i_nonz = ~(pheno_.isnan())

                self._lmm.setK(self.K[i_nonz,i_nonz])
                self._lmm.setSNPs(self.snps[i_nonz])
                self._lmm.setPheno(pheno_[i_nonz,np.newaxis])
                self._lmm.setCovs(self.covs[i_nonz])
                self._lmm.process()
                self.pvalues[ip:ip+1] = self._lmm.getPv()
                self.beta_snp[ip:ip+1] = self._lmm.getBetaSNP()
                self.beta_ste[ip:ip+1] = self._lmm.getBetaSNPste()
                self.ldelta_0[ip:ip+1] = self._lmm.getLdelta0()
                self.ldelta_alt[ip:ip+1] = self._lmm.getLdeltaAlt()
                self.NLL_0[ip:ip+1] = self._lmm.getNLL0()
                self.NLL_alt[ip:ip+1] = self._lmm.getNLLAlt()
                self.test_statistics[ip:ip+1] = self._lmm.getTestStatistics()
                pass
        if self._lmm.getTestStatistics() == self._lmm.TEST_LRT and self.test != "lrt":
            raise NotImplementedError("only f and lrt are implemented")
        elif self._lmm.getTestStatistics() == self._lmm.TEST_F and self.test != "f":
            raise NotImplementedError("only f and lrt are implemented")

        if self._lmm.getTestStatistics() == self._lmm.TEST_F:
            self.test_statistics = (self.beta_snp*self.beta_snp)/(self.beta_ste*self.beta_ste)
        if self._lmm.getTestStatistics() == self._lmm.TEST_LRT:
            self.test_statistics = 2.0 * (self.NLL_0 - self.NLL_alt)
        t1=time.time()

        if self.verbose:
            print(("finished GWAS testing in %.2f seconds" %(t1-t0)))

    def setCovs(self,covs):
        self._lmm.setCovs(covs)

    def getBetaSNP(self):
        """
        Returns:
            ndarray: [P, S] ndarray of SNP effect sizes.
        """
        return self.beta_snp

    def getPv(self):
        """
        Returns:
            ndarray: [P, S] ndarray of P-values.
        """
        return self.pvalues

    def getBetaSNPste(self):
        """
        Returns:
            ndarray: [P, S] ndarray of standard errors over SNP effects.
        """
        beta = self.getBetaSNP()
        pv = self.getPv()
        z = sp.sign(beta)*sp.sqrt(st.chi2(1).isf(pv))
        ste = beta/z
        return ste 

