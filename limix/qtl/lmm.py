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

import time

import numpy as np
import scipy as sp
import scipy.stats as st


class LMM:
    r"""Basic class for univariate single-variant association testing with LMMs.

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

    Examples
    --------

        Example of single-variant association testing using a linear mixed model.

        .. doctest::

            >>> from numpy.random import RandomState
            >>> from limix.qtl import LMM
            >>> from numpy import dot
            >>> from numpy import set_printoptions
            >>> set_printoptions(4)
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
            >>> lmm = LMM(snps, pheno, kinship)
            >>> pv = lmm.getPv()
            >>> beta = lmm.getBetaSNP()
            >>> beta_ste = lmm.getBetaSNPste()
            >>>
            >>> print(pv.shape)
            (1, 1000)
            >>> print(pv[:,:4])
            [[0.8571 0.4668 0.5872 0.5589]]
            >>> print(beta[:, :4])
            [[ 0.0436 -0.1695 -0.12    0.1388]]
            >>> print(beta_ste[:, :4])
            [[0.2422 0.2329 0.221  0.2375]]

        As shown in the exmaple below,
        the method can be applied directly on multiple phenotypes.
        Each phenotype is indipendently tested against all variants one-by-one.

        .. doctest::

            >>> random = RandomState(1)
            >>>
            >>> P = 3
            >>> phenos = random.randn(N, P)
            >>>
            >>> lmm = LMM(snps, phenos, kinship)
            >>> pv = lmm.getPv()
            >>>
            >>> print(pv.shape)
            (3, 1000)
            >>> print(pv[:,:4])
            [[0.4712 0.694  0.3369 0.5221]
             [0.7336 0.8799 0.8109 0.1071]
             [0.0662 0.9203 0.2873 0.8268]]
    """

    def __init__(
        self,
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
        # create column of 1 for fixed if nothing provide
        import limix_legacy.deprecated
        import limix_legacy.deprecated as dlimix_legacy

        if len(pheno.shape) == 1:
            pheno = pheno[:, sp.newaxis]

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
        self.N = self.pheno.shape[0]
        self.P = self.pheno.shape[1]
        self.Iok = ~(np.isnan(self.pheno).any(axis=1))
        if self.K is None:
            self.searchDelta = False
            self.K = np.eye(self.snps.shape[0])
        if self.covs is None:
            self.covs = np.ones((self.snps.shape[0], 1))

        self._lmm = None
        # run
        self.verbose = verbose
        self.process()

    def process(self):

        import limix_legacy.deprecated
        import limix_legacy.deprecated as dlimix_legacy

        t0 = time.time()
        if self._lmm is None:
            self._lmm = limix_legacy.deprecated.CLMM()
            self._lmm.setK(self.K)
            self._lmm.setSNPs(self.snps)
            self._lmm.setPheno(self.pheno)
            self._lmm.setCovs(self.covs)
            if self.test == "lrt":
                self._lmm.setTestStatistics(self._lmm.TEST_LRT)
            elif self.test == "f":
                self._lmm.setTestStatistics(self._lmm.TEST_F)
            else:
                print((self.test))
                raise NotImplementedError("only f and lrt are implemented")
            # set number of delta grid optimizations?
            self._lmm.setNumIntervals0(self.NumIntervalsDelta0)
            if self.searchDelta:
                self._lmm.setNumIntervalsAlt(self.NumIntervalsDeltaAlt)
            else:
                self._lmm.setNumIntervalsAlt(0)

        if not np.isnan(self.pheno).any():
            # process
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
                raise Exception(
                    "cannot reuse a CLMM object if missing variables are present"
                )
            else:
                self._lmm = limix_legacy.deprecated.CLMM()
            # test all phenotypes separately
            self.pvalues = np.zeros((self.phenos.shape[1], self.snps.shape[1]))
            self.beta_snp = np.zeros((self.phenos.shape[1], self.snps.shape[1]))
            self.beta_ste = np.zeros((self.phenos.shape[1], self.snps.shape[1]))
            self.ldelta_0 = np.zeros((self.phenos.shape[1], self.snps.shape[1]))
            self.ldelta_alt = np.zeros((self.phenos.shape[1], self.snps.shape[1]))
            self.NLL_0 = np.zeros((self.phenos.shape[1], self.snps.shape[1]))
            self.NLL_alt = np.zeros((self.phenos.shape[1], self.snps.shape[1]))
            self.test_statistics = np.zeros((self.phenos.shape[1], self.snps.shape[1]))
            for ip in np.arange(self.phenos.shape[1]):
                pheno_ = self.phenos[:, ip]
                i_nonz = ~(pheno_.isnan())

                self._lmm.setK(self.K[i_nonz, i_nonz])
                self._lmm.setSNPs(self.snps[i_nonz])
                self._lmm.setPheno(pheno_[i_nonz, np.newaxis])
                self._lmm.setCovs(self.covs[i_nonz])
                self._lmm.process()
                self.pvalues[ip : ip + 1] = self._lmm.getPv()
                self.beta_snp[ip : ip + 1] = self._lmm.getBetaSNP()
                self.beta_ste[ip : ip + 1] = self._lmm.getBetaSNPste()
                self.ldelta_0[ip : ip + 1] = self._lmm.getLdelta0()
                self.ldelta_alt[ip : ip + 1] = self._lmm.getLdeltaAlt()
                self.NLL_0[ip : ip + 1] = self._lmm.getNLL0()
                self.NLL_alt[ip : ip + 1] = self._lmm.getNLLAlt()
                self.test_statistics[ip : ip + 1] = self._lmm.getTestStatistics()
                pass
        if self._lmm.getTestStatistics() == self._lmm.TEST_LRT and self.test != "lrt":
            raise NotImplementedError("only f and lrt are implemented")
        elif self._lmm.getTestStatistics() == self._lmm.TEST_F and self.test != "f":
            raise NotImplementedError("only f and lrt are implemented")

        if self._lmm.getTestStatistics() == self._lmm.TEST_F:
            self.test_statistics = (self.beta_snp * self.beta_snp) / (
                self.beta_ste * self.beta_ste
            )
        if self._lmm.getTestStatistics() == self._lmm.TEST_LRT:
            self.test_statistics = 2.0 * (self.NLL_0 - self.NLL_alt)
        t1 = time.time()

        if self.verbose:
            print(("finished GWAS testing in %.2f seconds" % (t1 - t0)))

    def setCovs(self, covs):
        self._lmm.setCovs(covs)

    def getBetaSNP(self):
        """
        Returns:
            ndarray: (`P`, `S`) ndarray of SNP effect sizes.
        """
        return self.beta_snp

    def getPv(self):
        """
        Returns:
            ndarray: (`P`, `S`) ndarray of P-values.
        """
        return self.pvalues

    def getBetaSNPste(self):
        """
        Returns:
            ndarray: (`P`, `S`) ndarray of standard errors over SNP effects.
        """
        beta = self.getBetaSNP()
        pv = self.getPv()
        z = sp.sign(beta) * sp.sqrt(st.chi2(1).isf(pv))
        ste = beta / z
        return ste
