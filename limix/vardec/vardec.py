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

import copy
import sys
import warnings

import scipy as sp
import scipy.linalg
import scipy.stats


class VarianceDecomposition(object):
    r"""Class for variance decomposition

    Args:
        Y (ndarray):
            (`N`, `P`) phenotype matrix
        standardize (bool):
            if True, impute missing phenotype values by mean value,
            zero-mean and unit-variance phenotype.
            The default value is False.

    Examples
    --------

        Basic example of usage.

        .. doctest::

            >>> from numpy.random import RandomState
            >>> from limix.vardec import VarianceDecomposition
            >>> from numpy import dot, eye, ones
            >>> from numpy import set_printoptions
            >>> set_printoptions(4)
            >>> random = RandomState(1)
            >>>
            >>> N = 100
            >>> P = 3
            >>>
            >>> # generate data
            >>> pheno = random.randn(N, P)
            >>> W = random.randn(N, 10)
            >>> kinship = dot(W, W.T) / float(10)
            >>> kinship+= 1e-4*eye(N)
            >>> F1 = ones((N, 1))
            >>> A1 = eye(P)
            >>> F2 = random.randn(N, 2)
            >>> A2 = ones((1, P))
            >>>
            >>> vc = VarianceDecomposition(pheno)
            >>> vc.addFixedEffect(F1, A1)
            >>> vc.addFixedEffect(F2, A2)
            >>> vc.addRandomEffect(K=kinship)
            >>> vc.addRandomEffect(is_noise=True)
            >>> conv = vc.optimize()
            >>> print(vc.getTraitCovar(0)) # doctest: +SKIP
            [[ 0.0097  0.023  -0.0032]
             [ 0.023   0.0562 -0.0062]
             [-0.0032 -0.0062  0.0096]]
            >>> print(vc.getTraitCovar(1)) # doctest: +SKIP
            [[ 1.0029 -0.1212 -0.0052]
             [-0.1212  0.8284 -0.0412]
             [-0.0052 -0.0412  0.8147]]
            >>> print(vc.getWeights(0))
            [[0.039  0.0899 0.1213]]
            >>> print(vc.getWeights(1))
            [[-0.0176]
             [-0.0095]]
    """

    def __init__(self, Y, standardize=False):

        from limix.util import preprocess

        # check whether Y is a vector, if yes reshape
        if (len(Y.shape) == 1):
            Y = Y[:, sp.newaxis]

        # create column of 1 for fixed if nothing provided
        self.N = Y.shape[0]
        self.P = Y.shape[1]
        self.Nt = self.N * self.P
        self.Iok = ~(sp.isnan(Y).any(axis=1))

        # outsourced to handle missing values:
        if standardize:
            Y = preprocess.standardize(Y)

        # set properties
        self.Y = Y
        self.n_randEffs = 0
        self.noisPos = None

        # fixed and ranodm effects
        self.sample_covars = []
        self.sample_cross_covars = []
        self.trait_covars = []
        self.sample_designs = []
        self.sample_test_designs = []
        self.trait_designs = []

        # desinc
        self._desync()

        # TODO: check if necessary
        self.Ntest = None

    def setY(self, Y, standardize=False):
        """
        Set phenotype matrix

        Args:
            Y:              phenotype matrix [N, P]
            standardize:    if True, phenotype is standardized (zero mean, unit variance)
        """
        from limix.util import preprocess

        assert Y.shape[
            0] == self.N, 'VarianceDecomposition:: Incompatible shape'
        assert Y.shape[
            1] == self.P, 'VarianceDecomposition:: Incompatible shape'

        if standardize:
            Y = preprocess.standardize(Y)

        # check that missing values match the current structure
        assert (~(sp.isnan(Y).any(axis=1)) == self.Iok).all(
        ), 'VarianceDecomposition:: pattern of missing values needs to match Y given at initialization'

        # TODO: update the current code

        self.Y = Y
        self.vd.setPheno(Y)

        self.init = False
        self.optimum = None

        self.cache['Sigma'] = None
        self.cache['Hessian'] = None

    def setTestSampleSize(self, Ntest):
        """
        set sample size of the test dataset

        Args:
            Ntest:        sample size of the test set
        """
        self.Ntest = Ntest

    def addRandomEffect(self,
                        K=None,
                        is_noise=False,
                        normalize=False,
                        Kcross=None,
                        trait_covar_type='freeform',
                        rank=1,
                        fixed_trait_covar=None,
                        jitter=1e-4):
        """
        Add random effects term.

        Args:
            K:      Sample Covariance Matrix [N, N]
            is_noise:   Boolean indicator specifying if the matrix is homoscedastic noise (weighted identity covariance) (default False)
            normalize:  Boolean indicator specifying if K has to be normalized such that K.trace()=N.
            Kcross:            NxNtest cross covariance for predictions
            trait_covar_type: type of covaraince to use. Default 'freeform'. possible values are
                            'freeform':  general semi-definite positive matrix,
                            'fixed': fixed matrix specified in fixed_trait_covar,
                            'diag': diagonal matrix,
                            'lowrank': low rank matrix. The rank of the lowrank part is specified in the variable rank,
                            'lowrank_id': sum of a lowrank matrix and a multiple of the identity. The rank of the lowrank part is specified in the variable rank,
                            'lowrank_diag': sum of a lowrank and a diagonal matrix. The rank of the lowrank part is specified in the variable rank,
                            'block': multiple of a matrix of ones,
                            'block_id': sum of a multiple of a matrix of ones and a multiple of the idenity,
                            'block_diag': sum of a multiple of a matrix of ones and a diagonal matrix,
            rank:       rank of a possible lowrank component (default 1)
            fixed_trait_covar:   PxP matrix for the (predefined) trait-to-trait covariance matrix if fixed type is used
            jitter:        diagonal contribution added to trait-to-trait covariance matrices for regularization
        """
        assert K is not None or is_noise, 'VarianceDecomposition:: Specify covariance structure'

        if is_noise:
            assert self.noisPos is None, 'VarianceDecomposition:: Noise term already exists'
            K = sp.eye(self.N)
            Kcross = None
            self.noisPos = self.n_randEffs
        else:
            assert K.shape[
                0] == self.N, 'VarianceDecomposition:: Incompatible shape for K'
            assert K.shape[
                1] == self.N, 'VarianceDecomposition:: Incompatible shape for K'
        if Kcross is not None:
            assert self.Ntest is not None, 'VarianceDecomposition:: specify Ntest for predictions (method VarianceDecomposition::setTestSampleSize)'
            assert Kcross.shape[
                0] == self.N, 'VarianceDecomposition:: Incompatible shape for Kcross'
            assert Kcross.shape[
                1] == self.Ntest, 'VarianceDecomposition:: Incompatible shape for Kcross'

        from limix.util.preprocess import covar_rescaling_factor
        if normalize:
            cc = covar_rescaling_factor(K)
            K *= cc
            if Kcross is not None:
                Kcross *= cc

        # add random effect
        self.sample_covars.append(K)
        self.sample_cross_covars.append(Kcross)
        _cov = self._buildTraitCovar(
            trait_covar_type=trait_covar_type,
            rank=rank,
            fixed_trait_covar=fixed_trait_covar,
            jitter=jitter)
        self.trait_covars.append(_cov)
        self.n_randEffs += 1

        # TODO: define _sync function
        self._desync()

    def addFixedEffect(self, F=None, A=None, Ftest=None):
        """
        add fixed effect term to the model

        Args:
            F:     sample design matrix for the fixed effect [N,K]
            A:     trait design matrix for the fixed effect (e.g. sp.ones((1,P)) common effect; sp.eye(P) any effect) [L,P]
            Ftest: sample design matrix for test samples [Ntest,K]
        """
        if A is None:
            A = sp.eye(self.P)
        if F is None:
            F = sp.ones((self.N, 1))
            if self.Ntest is not None:
                Ftest = sp.ones((self.Ntest, 1))

        assert A.shape[
            1] == self.P, 'VarianceDecomposition:: A has incompatible shape'
        assert F.shape[
            0] == self.N, 'VarianceDecimposition:: F has incompatible shape'

        if Ftest is not None:
            assert self.Ntest is not None, 'VarianceDecomposition:: specify Ntest for predictions (method VarianceDecomposition::setTestSampleSize)'
            assert Ftest.shape[
                0] == self.Ntest, 'VarianceDecimposition:: Ftest has incompatible shape'
            assert Ftest.shape[1] == F.shape[
                1], 'VarianceDecimposition:: Ftest has incompatible shape'

        # add fixed effect
        self.sample_designs.append(F)
        self.sample_test_designs.append(Ftest)
        self.trait_designs.append(A)

        self._desync()

    def optimize(self,
                 init_method='default',
                 inference=None,
                 n_times=10,
                 perturb=False,
                 pertSize=1e-3,
                 verbose=False):
        """
        Train the model using the specified initialization strategy

        Args:
            init_method:    initialization strategy:
                                'default': variance is equally split across the different random effect terms. For mulit-trait models the empirical covariance between traits is used
                                'random': variance component parameters (scales) are sampled from a normal distribution with mean 0 and std 1,
                                None: no initialization is considered. Initial parameters can be specifies by using the single covariance getTraitCovarfun()
            inference:      inference gp method, by default algebrically efficient inference (i.e., gp2kronSum, gp2KronSumLR, gp3KronSumLR) will be used when possible.
                            For models with high a standard inference scheme (gp_base) will be used.
            n_times:        number of restarts to converge
            perturb:        if true, the initial point (if random initializaiton is not being used) is perturbed with gaussian noise for each restart (default, False)
            perturbSize:    std of the gassian noise used to perturb the initial point
            verbose:        print if convergence is achieved and how many restarted were needed
        """
        # verbose = dlimix.getVerbose(verbose)
        assert init_method in [
            'default', 'random', None
        ], 'VarianceDecomposition: specified init_method not valid'
        if init_method == 'default':
            self._init_params_default()

        # determine if repeats should be considered
        if init_method is not 'random' and not perturb:
            n_times = 1

        # determine inference
        if inference is None:
            inference = self._det_inference()
        else:
            self._check_inference(inference)
        self._inference = inference

        if self.gp is None:
            self._initGP()

        params0 = self.gp.getParams()
        for i in range(n_times):
            if init_method == 'random':
                params = {'covar': sp.randn(params0['covar'].shape[0])}
                self.gp.setParams(params)
            elif perturb:
                params = {
                    'covar':
                    params0['covar'] +
                    pertSize * sp.randn(params0['covar'].shape[0])
                }
                self.gp.setParams(params)
            conv, info = self.gp.optimize(verbose=verbose)
            if conv:
                break

        if verbose:
            if not conv:
                print(
                    'No local minimum found for the tested initialization points'
                )
            else:
                print(('Local minimum found at iteration %d' % i))

        return conv

    def getLML(self):
        """
        Return log marginal likelihood
        """
        assert self.init, 'VarianceDecomposition:: GP not initialised'
        return self.gp.LML()

    def getLMLgrad(self):
        """
        Return gradient of log-marginal likelihood
        """
        assert self.init, 'VarianceDecomposition:: GP not initialised'
        return self.gp.LML_grad()['covar']

    def getWeights(self, term_i=None):
        """
        Return weights for fixed effect term term_i

        Args:
            term_i:     fixed effect term index
        Returns:
            weights of the spefied fixed effect term.
            The output will be a KxL matrix of weights will be returned,
            where K is F.shape[1] and L is A.shape[1] of the correspoding fixed effect term
            (L will be always 1 for single-trait analysis).
        """
        assert self.init, 'GP not initialised'
        if term_i is None:
            if self.gp.mean.n_terms == 1:
                term_i = 0
            else:
                print('VarianceDecomposition: Specify fixed effect term index')
        return self.gp.mean.B[term_i]

    def getTraitCovarFun(self, term_i):
        """
        Return the trait limix covariance function for term term_i

        Args:
            term_i:     index of the ranodm effect to consider
        Returns:
            LIMIX cvoariance function
        """
        assert term_i < self.n_randEffs, 'VarianceDecomposition:: specied term out of range'
        return self.trait_covars[term_i]

    def getTraitCovar(self, term_i=None):
        """
        Return the estimated trait covariance matrix for term_i (or the total if term_i is None)
            To retrieve the matrix of correlation coefficient use \see getTraitCorrCoef

        Args:
            term_i:     index of the random effect term we want to retrieve the covariance matrix
        Returns:
            estimated trait covariance
        """
        assert term_i < self.n_randEffs, 'VarianceDecomposition:: specied term out of range'

        if term_i is None:
            RV = sp.zeros((self.P, self.P))
            for term_i in range(self.n_randEffs):
                RV += self.getTraitCovarFun().K()
        else:
            assert term_i < self.n_randEffs, 'Term index non valid'
            RV = self.getTraitCovarFun(term_i).K()
        return RV

    def getTraitCorrCoef(self, term_i=None):
        """
        Return the estimated trait correlation coefficient matrix for term_i (or the total if term_i is None)
            To retrieve the trait covariance matrix use \see getTraitCovar

        Args:
            term_i:     index of the random effect term we want to retrieve the correlation coefficients
        Returns:
            estimated trait correlation coefficient matrix
        """
        cov = self.getTraitCovar(term_i)
        stds = sp.sqrt(cov.diagonal())[:, sp.newaxis]
        RV = cov / stds / stds.T
        return RV

    def getVarianceComps(self, univariance=False):
        """
        Return the estimated variance components

        Args:
            univariance:   Boolean indicator, if True variance components are normalized to sum up to 1 for each trait
        Returns:
            variance components of all random effects on all phenotypes [P, n_randEffs matrix]
        """
        RV = sp.zeros((self.P, self.n_randEffs))
        for term_i in range(self.n_randEffs):
            RV[:, term_i] = self.getTraitCovar(term_i).diagonal()
        if univariance:
            RV /= RV.sum(1)[:, sp.newaxis]
        return RV

    ##################################
    # INTERNAL METHODS
    #################################
    def _desync(self):
        self._inference = None
        self.gp = None

    @property
    def init(self):
        return self.gp is not None

    def _init_params_default(self):
        """
        Internal method for default parameter initialization
        """
        # if there are some nan -> mean impute
        Yimp = self.Y.copy()
        Inan = sp.isnan(Yimp)
        Yimp[Inan] = Yimp[~Inan].mean()
        if self.P == 1:
            C = sp.array([[Yimp.var()]])
        else:
            C = sp.cov(Yimp.T)
        C /= float(self.n_randEffs)
        for ti in range(self.n_randEffs):
            self.getTraitCovarFun(ti).setCovariance(C)

    def _det_inference(self):
        """
        Internal method for determining the inference method
        """
        # 2 random effects with complete design -> gp2KronSum
        # TODO: add check for low-rankness, use GP3KronSumLR and GP2KronSumLR
        # when possible
        if (self.n_randEffs == 2) and (~sp.isnan(self.Y).any()):
            rv = 'GP2KronSum'
        else:
            rv = 'GP'
        return rv

    def _check_inference(self, inference):
        """
        Internal method for checking that the selected inference scheme is compatible with the specified model
        """
        if inference == 'GP2KronSum':
            assert self.n_randEffs == 2, 'VarianceDecomposition: for fast inference number of random effect terms must be == 2'
            assert not sp.isnan(self.Y).any(
            ), 'VarianceDecomposition: fast inference available only for complete phenotype designs'
        # TODO: add GP3KronSumLR, GP2KronSumLR

    def _initGP(self):
        """
        Internal method for initialization of the GP inference objetct
        """
        from limix.util.util_functions import vec
        from limix_core.mean import MeanKronSum
        from limix_core.gp import GP2KronSum
        from limix_core.covar import SumCov
        if self._inference == 'GP2KronSum':
            signalPos = sp.where(
                sp.arange(self.n_randEffs) != self.noisPos)[0][0]
            gp = GP2KronSum(
                Y=self.Y,
                F=self.sample_designs,
                A=self.trait_designs,
                Cg=self.trait_covars[signalPos],
                Cn=self.trait_covars[self.noisPos],
                R=self.sample_covars[signalPos])
        else:
            mean = MeanKronSum(self.Y, self.sample_designs, self.trait_designs)
            Iok = vec(~sp.isnan(mean.Y))[:, 0]
            if Iok.all():
                Iok = None
            covar = SumCov(*[
                KronCov(self.trait_covars[i], self.sample_covars[i], Iok=Iok)
                for i in range(self.n_randEffs)
            ])
            gp = GP(covar=covar, mean=mean)
        self.gp = gp

    def _buildTraitCovar(self,
                         trait_covar_type='freeform',
                         rank=1,
                         fixed_trait_covar=None,
                         jitter=1e-4):
        """
        Internal functions that builds the trait covariance matrix using the LIMIX framework

        Args:
            trait_covar_type: type of covaraince to use. Default 'freeform'. possible values are
            rank:                rank of a possible lowrank component (default 1)
            fixed_trait_covar:   PxP matrix for the (predefined) trait-to-trait covariance matrix if fixed type is used
            jitter:              diagonal contribution added to freeform covariance matrices for regularization
        Returns:
            LIMIX::Covariance for Trait covariance matrix
        """
        from limix_core.covar import (FreeFormCov, FixedCov, DiagonalCov,
                                      LowRankCov, SumCov)

        assert trait_covar_type in [
            'freeform', 'diag', 'lowrank', 'lowrank_id', 'lowrank_diag',
            'block', 'block_id', 'block_diag', 'fixed'
        ], 'VarianceDecomposition:: trait_covar_type not valid'

        if trait_covar_type == 'freeform':
            cov = FreeFormCov(self.P, jitter=jitter)
        elif trait_covar_type == 'fixed':
            assert fixed_trait_covar is not None, 'VarianceDecomposition:: set fixed_trait_covar'
            assert fixed_trait_covar.shape[
                0] == self.P, 'VarianceDecomposition:: Incompatible shape for fixed_trait_covar'
            assert fixed_trait_covar.shape[
                1] == self.P, 'VarianceDecomposition:: Incompatible shape for fixed_trait_covar'
            cov = FixedCov(fixed_trait_covar)
        elif trait_covar_type == 'diag':
            cov = DiagonalCov(self.P)
        elif trait_covar_type == 'lowrank':
            cov = LowRankCov(self.P, rank=rank)
        elif trait_covar_type == 'lowrank_id':
            cov = SumCov(
                LowRankCov(self.P, rank=rank), FixedCov(sp.eye(self.P)))
        elif trait_covar_type == 'lowrank_diag':
            cov = SumCov(LowRankCov(self.P, rank=rank), DiagonalCov(self.P))
        elif trait_covar_type == 'block':
            cov = FixedCov(sp.ones([self.P, self.P]))
        elif trait_covar_type == 'block_id':
            cov1 = FixedCov(sp.ones([self.P, self.P]))
            cov2 = FixedCov(sp.eye(self.P))
            cov = SumCov(cov1, cov2)
        elif trait_covar_type == 'block_diag':
            cov1 = FixedCov(sp.ones([self.P, self.P]))
            cov2 = FixedCov(sp.eye(self.P))
            cov = SumCov(cov1, cov2)
        return cov

    #####################################
    # METHODS TO UPDATE (OR DEPRECATED)
    #####################################

    def optimize_with_repeates(self,
                               fast=None,
                               verbose=None,
                               n_times=10,
                               lambd=None,
                               lambd_g=None,
                               lambd_n=None):
        """
        Train the model repeadly up to a number specified by the users with random restarts and
        return a list of all relative minima that have been found. This list is sorted according to
        least likelihood. Each list term is a dictionary with keys "counter", "LML", and "scales".

        After running this function, the vc object will be set at the last iteration. Thus, if you
        wish to get the vc object of one of the repeats, then set the scales. For example:

        vc.setScales(scales=optimize_with_repeates_output[0]["scales"])

        Args:
            fast:       Boolean. if set to True initalize kronSumGP
            verbose:    Boolean. If set to True, verbose output is produced. (default True)
            n_times:    number of re-starts of the optimization. (default 10)
        """
        verbose = dlimix.getVerbose(verbose)

        if not self.init:
            self._initGP(fast)

        opt_list = []

        fixed0 = sp.zeros_like(self.gp.getParams()['dataTerm'])

        # minimize n_times
        for i in range(n_times):

            scales1 = self._getScalesRand()
            fixed1 = 1e-1 * sp.randn(fixed0.shape[0], fixed0.shape[1])
            conv = self.trainGP(
                fast=fast,
                scales0=scales1,
                fixed0=fixed1,
                lambd=lambd,
                lambd_g=lambd_g,
                lambd_n=lambd_n)

            if conv:
                # compare with previous minima
                temp = 1
                for j in range(len(opt_list)):
                    if sp.allclose(
                            abs(self.getScales()), abs(opt_list[j]['scales'])):
                        temp = 0
                        opt_list[j]['counter'] += 1
                        break
                if temp == 1:
                    opt = {}
                    opt['counter'] = 1
                    opt['LML'] = self.getLML()
                    opt['scales'] = self.getScales()
                    opt_list.append(opt)

        # sort by LML
        LML = sp.array([opt_list[i]['LML'] for i in range(len(opt_list))])
        index = LML.argsort()[::-1]
        out = []
        if verbose:
            print("\nLocal mimima\n")
            print("n_times\t\tLML")
            print("------------------------------------")

        for i in range(len(opt_list)):
            out.append(opt_list[index[i]])
            if verbose:
                print(("%d\t\t%f" % (opt_list[index[i]]['counter'],
                                     opt_list[index[i]]['LML'])))
                print("")

        return out

    """ standard errors """

    def getTraitCovarStdErrors(self, term_i):
        """
        Returns standard errors on trait covariances from term_i (for the covariance estimate \see getTraitCovar)

        Args:
            term_i:     index of the term we are interested in
        """
        assert self.init, 'GP not initialised'
        assert self.fast is False, 'Not supported for fast implementation'

        if self.P == 1:
            out = (2 * self.getScales()[term_i])**2 * \
                self._getLaplaceCovar()[term_i, term_i]
        else:
            C = self.vd.getTerm(term_i).getTraitCovar()
            n_params = C.getNumberParams()
            par_index = 0
            for term in range(term_i - 1):
                par_index += self.vd.getTerm(term_i).getNumberScales()
            Sigma1 = self._getLaplaceCovar()[par_index:(
                par_index + n_params), :][:, par_index:(par_index + n_params)]
            out = sp.zeros((self.P, self.P))
            for param_i in range(n_params):
                out += C.Kgrad_param(param_i)**2 * Sigma1[param_i, param_i]
                for param_j in range(param_i):
                    out += 2 * abs(
                        C.Kgrad_param(param_i) * C.Kgrad_param(param_j)
                    ) * Sigma1[param_i, param_j]
        out = sp.sqrt(out)
        return out

    def getVarianceCompStdErrors(self, univariance=False):
        """
        Return the standard errors on the estimated variance components (for variance component estimates \see getVarianceComps)

        Args:
            univariance:   Boolean indicator, if True variance components are normalized to sum up to 1 for each trait
        Returns:
            standard errors on variance components [P, n_randEffs matrix]
        """
        RV = sp.zeros((self.P, self.n_randEffs))
        for term_i in range(self.n_randEffs):
            # RV[:,term_i] = self.getTraitCovarStdErrors(term_i).diagonal()
            RV[:, term_i] = self.getTraitCovarStdErrors(term_i)
        var = self.getVarianceComps()
        if univariance:
            RV /= var.sum(1)[:, sp.newaxis]
        return RV

    """
    CODE FOR PREDICTIONS
    """

    def predictPhenos(self, use_fixed=None, use_random=None):
        """
        predict the conditional mean (BLUP)

        Args:
            use_fixed:        list of fixed effect indeces to use for predictions
            use_random:        list of random effect indeces to use for predictions
        Returns:
            predictions (BLUP)
        """
        assert self.noisPos is not None, 'No noise element'
        assert self.init, 'GP not initialised'
        assert self.Ntest is not None, 'VarianceDecomposition:: specify Ntest for predictions (method VarianceDecomposition::setTestSampleSize)'

        use_fixed = list(range(self.n_fixedEffs))
        use_random = list(range(self.n_randEffs))

        KiY = self.gp.agetKEffInvYCache()

        if not self.fast:
            KiY = KiY.reshape(self.P, self.N).T

        Ypred = sp.zeros((self.Ntest, self.P))

        # predicting from random effects
        for term_i in use_random:
            if term_i != self.noisPos:
                Kstar = self.Kstar[term_i]
                if Kstar is None:
                    warnings.warn(
                        'warning: random effect term %d not used for predictions as it has None cross covariance'
                        % term_i)
                    continue
                term = sp.dot(Kstar.T, KiY)
                if self.P > 1:
                    C = self.getTraitCovar(term_i)
                    term = sp.dot(term, C)
                else:
                    term *= self.getVarianceComps()[0, term_i]
                Ypred += term

        # predicting from fixed effects
        weights = self.getWeights()
        w_i = 0
        for term_i in use_fixed:
            Fstar = self.Fstar[term_i]
            if Fstar is None:
                warnings.warn(
                    'warning: fixed effect term %d not used for predictions as it has None test sample design'
                    % term_i)
                continue
            if self.P == 1:
                A = sp.eye(1)
            else:
                A = self.vd.getDesign(term_i)
            Fstar = self.Fstar[term_i]
            W = weights[w_i:w_i + A.shape[0], 0:1].T
            term = sp.dot(Fstar, sp.dot(W, A))
            w_i += A.shape[0]
            Ypred += term

        return Ypred

    def crossValidation(self,
                        seed=0,
                        n_folds=10,
                        fullVector=True,
                        verbose=None,
                        D=None,
                        **keywords):
        """
        Split the dataset in n folds, predict each fold after training the model on all the others

        Args:
            seed:        seed
            n_folds:     number of folds to train the model on
            fullVector:  Bolean indicator, if true it stops if no convergence is observed for one of the folds, otherwise goes through and returns a pheno matrix with missing values
            verbose:     if true, prints the fold that is being used for predicitons
            **keywords:  params to pass to the function optimize
        Returns:
            Matrix of phenotype predictions [N,P]
        """
        verbose = dlimix.getVerbose(verbose)

        # split samples into training and test
        sp.random.seed(seed)
        r = sp.random.permutation(self.Y.shape[0])
        nfolds = 10
        Icv = sp.floor(((sp.ones(
            (self.Y.shape[0])) * nfolds) * r) / self.Y.shape[0])

        RV = {}
        if self.P == 1:
            RV['var'] = sp.zeros((nfolds, self.n_randEffs))
        else:
            RV['var'] = sp.zeros((nfolds, self.P, self.n_randEffs))

        Ystar = sp.zeros_like(self.Y)

        for fold_j in range(n_folds):

            if verbose:
                print((".. predict fold %d" % fold_j))

            Itrain = Icv != fold_j
            Itest = Icv == fold_j
            Ytrain = self.Y[Itrain, :]
            Ytest = self.Y[Itest, :]
            vc = VarianceDecomposition(Ytrain)
            vc.setTestSampleSize(Itest.sum())
            for term_i in range(self.n_fixedEffs):
                F = self.vd.getFixed(term_i)
                Ftest = F[Itest, :]
                Ftrain = F[Itrain, :]
                if self.P > 1:
                    A = self.vd.getDesign(term_i)
                else:
                    A = None
                vc.addFixedEffect(F=Ftrain, Ftest=Ftest, A=A)
            for term_i in range(self.n_randEffs):
                if self.P > 1:
                    tct = self.trait_covar_type[term_i]
                    rank = self.rank[term_i]
                    ftc = self.fixed_tc[term_i]
                    jitt = self.jitter[term_i]
                    if tct == 'lowrank_diag1' or tct == 'freeform1':
                        d = D[fold_j, :, term_i]
                    else:
                        d = None
                else:
                    tct = None
                    rank = None
                    ftc = None
                    jitt = None
                    d = None
                if term_i == self.noisPos:
                    vc.addRandomEffect(
                        is_noise=True,
                        trait_covar_type=tct,
                        rank=rank,
                        jitter=jitt,
                        fixed_trait_covar=ftc,
                        d=d)
                else:
                    R = self.vd.getTerm(term_i).getK()
                    Rtrain = R[Itrain, :][:, Itrain]
                    Rcross = R[Itrain, :][:, Itest]
                    vc.addRandomEffect(
                        K=Rtrain,
                        Kcross=Rcross,
                        trait_covar_type=tct,
                        rank=rank,
                        jitter=jitt,
                        fixed_trait_covar=ftc,
                        d=d)
            conv = vc.optimize(verbose=False, **keywords)
            if self.P == 1:
                RV['var'][fold_j, :] = vc.getVarianceComps()[0, :]
            else:
                RV['var'][fold_j, :, :] = vc.getVarianceComps()

            if fullVector:
                assert conv, 'VarianceDecompositon:: not converged for fold %d. Stopped here' % fold_j
            if conv:
                Ystar[Itest, :] = vc.predictPhenos()
            else:
                warnings.warn('not converged for fold %d' % fold_j)
                Ystar[Itest, :] = sp.nan

        return Ystar, RV

    """ INTERNAL FUNCIONS FOR PARAMETER INITIALIZATION """

    def _getH2singleTrait(self, K, verbose=None):
        """
        Internal function for parameter initialization
        estimate variance components and fixed effect using a linear mixed model with an intercept and 2 random effects (one is noise)
        Args:
            K:        covariance matrix of the non-noise random effect term
        """
        verbose = dlimix.getVerbose(verbose)
        # Fit single trait model
        varg = sp.zeros(self.P)
        varn = sp.zeros(self.P)
        fixed = sp.zeros((1, self.P))

        for p in range(self.P):
            y = self.Y[:, p:p + 1]
            # check if some sull value
            I = sp.isnan(y[:, 0])
            if I.sum() > 0:
                y = y[~I, :]
                _K = K[~I, :][:, ~I]
            else:
                _K = copy.copy(K)
            lmm = dlimix.CLMM()
            lmm.setK(_K)
            lmm.setSNPs(sp.ones((y.shape[0], 1)))
            lmm.setPheno(y)
            lmm.setCovs(sp.zeros((y.shape[0], 1)))
            lmm.setVarcompApprox0(-20, 20, 1000)
            lmm.process()
            delta = sp.exp(lmm.getLdelta0()[0, 0])
            Vtot = sp.exp(lmm.getLSigma()[0, 0])

            varg[p] = Vtot
            varn[p] = delta * Vtot
            fixed[:, p] = lmm.getBetaSNP()

            if verbose:
                print(p)

        sth = {}
        sth['varg'] = varg
        sth['varn'] = varn
        sth['fixed'] = fixed

        return sth

    def _getScalesDiag(self, termx=0):
        """
        Internal function for parameter initialization
        Uses 2 term single trait model to get covar params for initialization

        Args:
            termx:      non-noise term terms that is used for initialization
        """
        assert self.P > 1, 'VarianceDecomposition:: diagonal init_method allowed only for multi trait models'
        assert self.noisPos is not None, 'VarianceDecomposition:: noise term has to be set'
        assert termx < self.n_randEffs - 1, 'VarianceDecomposition:: termx>=n_randEffs-1'
        assert self.trait_covar_type[self.noisPos] not in [
            'lowrank', 'block', 'fixed'
        ], 'VarianceDecomposition:: diagonal initializaiton not posible for such a parametrization'
        assert self.trait_covar_type[termx] not in [
            'lowrank', 'block', 'fixed'
        ], 'VarianceDecimposition:: diagonal initializaiton not posible for such a parametrization'
        scales = []
        res = self._getH2singleTrait(self.vd.getTerm(termx).getK())
        scaleg = sp.sqrt(res['varg'].mean())
        scalen = sp.sqrt(res['varn'].mean())
        for term_i in range(self.n_randEffs):
            if term_i == termx:
                _scales = scaleg * self.diag[term_i]
            elif term_i == self.noisPos:
                _scales = scalen * self.diag[term_i]
            else:
                _scales = 0. * self.diag[term_i]
            if self.jitter[term_i] > 0:
                _scales = sp.concatenate(
                    (_scales, sp.array([sp.sqrt(self.jitter[term_i])])))
            scales.append(_scales)
        return sp.concatenate(scales)

    def _getScalesPairwise(self, verbose=False, initDiagonal=False):
        """
        Internal function for parameter initialization
        Uses a single trait model for initializing variances and
        a pairwise model to initialize correlations
        """
        var = sp.zeros((self.P, 2))

        if initDiagonal:
            # 1. fit single trait model
            if verbose:
                print('.. fit single-trait model for initialization')
            vc = VarianceDecomposition(self.Y[:, 0:1])
            for term_i in range(self.n_randEffs):
                if term_i == self.noisPos:
                    vc.addRandomEffect(is_noise=True)
                else:
                    K = self.vd.getTerm(term_i).getK()
                    vc.addRandomEffect(K=K)
            scales0 = sp.sqrt(0.5) * sp.ones(2)

            for p in range(self.P):
                if verbose:
                    print(('   .. trait %d' % p))
                vc.setY(self.Y[:, p:p + 1])
                conv = vc.optimize(scales0=scales0)
                if not conv:
                    print('warning initialization not converged')
                var[p, :] = vc.getVarianceComps()[0, :]

        elif fastlmm_present:
            if verbose:
                print(
                    '.. fit single-trait model for initialization (using fastlmm)'
                )
            for p in range(self.P):
                if verbose:
                    print(('   .. trait %d' % p))
                covariates = None
                for term_i in range(self.n_randEffs):
                    if term_i == self.noisPos:
                        pass
                    else:
                        K = self.vd.getTerm(term_i).getK()
                varY = sp.var(self.Y[:, p:p + 1])
                lmm = fastLMM(X=covariates, Y=self.Y[:, p:p + 1], G=None, K=K)
                opt = lmm.findH2(nGridH2=100)
                h2 = opt['h2']
                var[p, :] = h2 * varY
                var[p, self.noisPos] = (1.0 - h2) * varY
        else:
            if verbose:
                print('.. random initialization of diagonal')
            var = sp.random.randn(var.shape[0], var.shape[1])
            var = var * var + 0.001
        # 2. fit pairwise model
        if verbose:
            print('.. fit pairwise model for initialization')
        vc = VarianceDecomposition(self.Y[:, 0:2])
        for term_i in range(self.n_randEffs):
            if term_i == self.noisPos:
                vc.addRandomEffect(is_noise=True, trait_covar_type='freeform')
            else:
                K = self.vd.getTerm(term_i).getK()
                vc.addRandomEffect(K=K, trait_covar_type='freeform')
        rho_g = sp.ones((self.P, self.P))
        rho_n = sp.ones((self.P, self.P))
        for p1 in range(self.P):
            for p2 in range(p1):
                if verbose:
                    print(('   .. fit pair (%d,%d)' % (p1, p2)))
                vc.setY(self.Y[:, [p1, p2]])
                scales0 = sp.sqrt(
                    sp.array([
                        var[p1, 0], 1e-4, var[p2, 0], 1e-4, var[p1, 1], 1e-4,
                        var[p2, 1], 1e-4
                    ]))
                conv = vc.optimize(scales0=scales0)
                if not conv:
                    print('warning initialization not converged')
                Cg = vc.getTraitCovar(0)
                Cn = vc.getTraitCovar(1)
                rho_g[p1, p2] = Cg[0, 1] / sp.sqrt(Cg.diagonal().prod())
                rho_n[p1, p2] = Cn[0, 1] / sp.sqrt(Cn.diagonal().prod())
                rho_g[p2, p1] = rho_g[p1, p2]
                rho_n[p2, p1] = rho_n[p1, p2]
        # 3. init
        Cg0 = rho_g * sp.dot(sp.sqrt(var[:, 0:1]), sp.sqrt(var[:, 0:1].T))
        Cn0 = rho_n * sp.dot(sp.sqrt(var[:, 1:2]), sp.sqrt(var[:, 1:2].T))
        offset_g = abs(sp.minimum(sp.linalg.eigh(Cg0)[0].min(), 0)) + 1e-4
        offset_n = abs(sp.minimum(sp.linalg.eigh(Cn0)[0].min(), 0)) + 1e-4
        Cg0 += offset_g * sp.eye(self.P)
        Cn0 += offset_n * sp.eye(self.P)
        Lg = sp.linalg.cholesky(Cg0)
        Ln = sp.linalg.cholesky(Cn0)
        Cg_params0 = sp.concatenate([Lg[:, p][:p + 1] for p in range(self.P)])
        Cn_params0 = sp.concatenate([Ln[:, p][:p + 1] for p in range(self.P)])
        scales0 = sp.concatenate(
            [Cg_params0, 1e-2 * sp.ones(1), Cn_params0, 1e-2 * sp.ones(1)])

        return scales0

    def _getScalesRand(self):
        """
        Internal function for parameter initialization
        Return a vector of random scales
        """
        if self.P > 1:
            scales = []
            for term_i in range(self.n_randEffs):
                _scales = sp.randn(self.diag[term_i].shape[0])
                if self.jitter[term_i] > 0:
                    _scales = sp.concatenate(
                        (_scales, sp.array([sp.sqrt(self.jitter[term_i])])))
                scales.append(_scales)
            scales = sp.concatenate(scales)
        else:
            scales = sp.randn(self.vd.getNumberScales())
        return scales

    def _perturbation(self):
        """
        Internal function for parameter initialization
        Returns Gaussian perturbation
        """
        if self.P > 1:
            scales = []
            for term_i in range(self.n_randEffs):
                _scales = sp.randn(self.diag[term_i].shape[0])
                if self.jitter[term_i] > 0:
                    _scales = sp.concatenate((_scales, sp.zeros(1)))
                scales.append(_scales)
            scales = sp.concatenate(scales)
        else:
            scales = sp.randn(self.vd.getNumberScales())
        return scales

    """ INTERNAL FUNCIONS FOR ESTIMATING PARAMETERS UNCERTAINTY """

    def _getHessian(self):
        """
        Internal function for estimating parameter uncertainty
        COMPUTES OF HESSIAN OF E(\theta) = - log L(\theta | X, y)
        """
        assert self.init, 'GP not initialised'
        assert self.fast is False, 'Not supported for fast implementation'

        if self.cache['Hessian'] is None:
            ParamMask = self.gp.getParamMask()['covar']
            std = sp.zeros(ParamMask.sum())
            H = self.gp.LMLhess_covar()
            It = (ParamMask[:, 0] == 1)
            self.cache['Hessian'] = H[It, :][:, It]

        return self.cache['Hessian']

    def _getLaplaceCovar(self):
        """
        Internal function for estimating parameter uncertainty
        Returns:
            the
        """
        assert self.init, 'GP not initialised'
        assert self.fast is False, 'Not supported for fast implementation'

        if self.cache['Sigma'] is None:
            self.cache['Sigma'] = sp.linalg.inv(self._getHessian())
        return self.cache['Sigma']

    def _getFisher(self):
        """
        Return the fisher information matrix over the parameters
        TO CHECK: IT DOES NOT WORK PROPERLY
        """
        Ctot = self.vd.getGP().getCovar()
        Ki = sp.linalg.inv(Ctot.K())
        n_scales = self.vd.getNumberScales()
        out = sp.zeros((n_scales, n_scales))
        for m in range(n_scales):
            out[m, m] = 0.5 * sp.trace(
                sp.dot(Ki,
                       sp.dot(
                           Ctot.Kgrad_param(m), sp.dot(Ki,
                                                       Ctot.Kgrad_param(m)))))
            for n in range(m):
                out[m, n] = 0.5 * sp.trace(
                    sp.dot(Ki,
                           sp.dot(
                               Ctot.Kgrad_param(m),
                               sp.dot(Ki, Ctot.Kgrad_param(n)))))
                out[n, m] = out[m, n]
        return out

    """ DEPRECATED """

    def _getModelPosterior(self, min):
        """
        USES LAPLACE APPROXIMATION TO CALCULATE THE BAYESIAN MODEL POSTERIOR
        """
        Sigma = self._getLaplaceCovar(min)
        n_params = self.vd.getNumberScales()
        ModCompl = 0.5 * n_params * \
            sp.log(2 * sp.pi) + 0.5 * sp.log(sp.linalg.det(Sigma))
        RV = min['LML'] + ModCompl
        return RV

    def trainGP(self,
                fast=None,
                scales0=None,
                fixed0=None,
                lambd=None,
                lambd_g=None,
                lambd_n=None):
        print('DEPRECATED')
