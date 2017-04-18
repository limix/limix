import sys
from limix.utils.preprocess import remove_dependent_cols
from limix.utils.util_functions import smartDumpDictHdf5

# core
from limix.core.gp import GP2KronSum
from limix.core.gp import GP2KronSumLR
from limix.core.gp import GP3KronSumLR
from limix.core.covar import FreeFormCov

import h5py
import pdb
import scipy as sp
import scipy.linalg as la
import scipy.stats as st
import time as TIME
import copy
import warnings
import os


class MTSet():
    """
    This class eases implementation of efficient set test algorithms.
    mtSet can account for relatedness and can be used for single- and multi-trait analysis
    """

    def __init__(self, Y=None, R=None, S_R=None, U_R=None, traitID=None, F=None, rank = 1):
        """
        Args:
            Y:          [N, P] phenotype matrix
            F:          [N, K] matrix of fixed effect design.
                        K is the number of per-trait covariates.
            R:          [N, N] genetic relatedness matrix between individuals.
                        In alternative to R, S_R and U_R can be provided.
                        If not specified a model without relatedness component is considered.
            S_R:        N vector of eigenvalues of R
            U_R:        [N, N] eigenvector matrix of R
            traiID:     P vector of the IDs of the phenotypes to analyze (optional)
            rank:       rank of the trait covariance matrix of the variance component to be tested (default is 1)
        """
        # data
        noneNone = S_R is not None and U_R is not None
        self.bgRE = R is not None or noneNone
        # fixed effect
        msg = 'The current implementation of the full rank mtSet'
        msg+= ' does not support covariates.'
        msg+= ' We reccommend to regress out covariates and'
        msg+= ' subsequently quantile normalize the phenotypes'
        msg+= ' to a normal distribution prior to use mtSet.'
        msg+= ' This can be done within the LIMIX framework using'
        msg+= ' the methods limix.utils.preprocess.regressOut and'
        msg+= ' limix.utils.preprocess.gaussianize'
        assert not (F is not None and self.bgRE), msg
        if F is not None:
            F = remove_dependent_cols(F)
            A = sp.eye(Y.shape[1])
        else:
            A = None
        #traitID
        if traitID is None:
            traitID = sp.array(['trait %d' % p for p in range(Y.shape[1])])
        self.setTraitID(traitID)
        #init covariance matrices and gp
        Cg = FreeFormCov(Y.shape[1])
        Cn = FreeFormCov(Y.shape[1])
        G  = 1. * (sp.rand(Y.shape[0],1)<0.2)
        if self.bgRE:
            self._gp = GP3KronSumLR(Y=Y, Cg=Cg, Cn=Cn, R=R, S_R=S_R, U_R=U_R, G=G, rank = rank)
        else:
            self._gp = GP2KronSumLR(Y=Y, Cn=Cn, G=G, F=F, A=A)
        # null model params
        self.null = None
        # calls itself for column-by-column trait analysis
        self.stSet = None
        self.nullST = None
        self.infoOpt   = None
        self.infoOptST = None
        pass

    ##################################################
    # Properties
    ##################################################
    @property
    def N(self):    return self._gp.covar.dim_r

    @property
    def P(self):    return self._gp.covar.dim_c

    @property
    def Cr(self):   return self._gp.covar.Cr

    @property
    def Cg(self):   return self._gp.covar.Cg

    @property
    def Cn(self):   return self._gp.covar.Cn

    @property
    def Y(self):    return self._gp.mean.Y

    @property
    def F(self):
        try:
            return self._gp.mean.F[0]
        except:
            return None

    @property
    def A(self):
        try:
            return self._gp.mean.A[0]
        except:
            return None

    @property
    def S_R(self):
        if self.bgRE:   RV = self._gp.covar.Sr()
        else:           RV = None
        return RV

    @property
    def U_R(self):
        if self.bgRE:   RV = self._gp.covar.Lr().T
        else:           RV = None
        return RV

    #########################################
    # Setters
    #########################################
    @Y.setter
    def Y(self, value):
        assert value.shape[0]==self.N, 'Dimension mismatch'
        assert value.shape[1]==self.P, 'Dimension mismatch'
        self._gp.mean.Y = Y

    def setTraitID(self,traitID):
        self.traitID = traitID

    ################################################
    # Fitting null model
    ###############################################
    def fitNull(self, verbose=False, cache=False, out_dir='./cache', fname=None, rewrite=False, seed=None, n_times=10, factr=1e3, init_method=None):
        """
        Fit null model
        """
        if seed is not None:    sp.random.seed(seed)

        read_from_file = False
        if cache:
            assert fname is not None, 'MultiTraitSetTest:: specify fname'
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            out_file = os.path.join(out_dir,fname)
            read_from_file = os.path.exists(out_file) and not rewrite

        RV = {}
        if read_from_file:
            f = h5py.File(out_file,'r')
            for key in list(f.keys()):
                RV[key] = f[key][:]
            f.close()
            self.setNull(RV)
        else:
            start = TIME.time()
            if self.bgRE:
                self._gpNull = GP2KronSum(Y=self.Y, F=None, A=None, Cg=self.Cg, Cn=self.Cn, R=None, S_R=self.S_R, U_R=self.U_R)
            else:
                self._gpNull = GP2KronSumLR(self.Y, self.Cn, G=sp.ones((self.N,1)), F=self.F, A=self.A)
                # freezes Cg to 0
                n_params = self._gpNull.covar.Cr.getNumberParams()
                self._gpNull.covar.Cr.setParams(1e-9 * sp.ones(n_params))
                self._gpNull.covar.act_Cr = False
            for i in range(n_times):
                params0 = self._initParams(init_method=init_method)
                self._gpNull.setParams(params0)
                conv, info = self._gpNull.optimize(verbose=verbose, factr=factr)
                if conv: break
            if not conv:    warnings.warn("not converged")
            LMLgrad = (self._gpNull.LML_grad()['covar']**2).mean()
            LML = self._gpNull.LML()
            if self._gpNull.mean.n_terms==1:
                RV['B'] = self._gpNull.mean.B[0]
            elif self._gpNull.mean.n_terms>1:
                warning.warn('generalize to more than 1 fixed effect term')
            if self.bgRE:
                RV['params0_g'] = self.Cg.getParams()
            else:
                RV['params0_g'] = sp.zeros_like(self.Cn.getParams())
            RV['params0_n'] = self.Cn.getParams()
            if self.bgRE:
                RV['Cg'] = self.Cg.K()
            else:
                RV['Cg'] = sp.zeros_like(self.Cn.K())
            RV['Cn'] = self.Cn.K()
            RV['conv'] = sp.array([conv])
            RV['time'] = sp.array([TIME.time()-start])
            RV['NLL0'] = sp.array([LML])
            RV['LMLgrad'] = sp.array([LMLgrad])
            RV['nit'] = sp.array([info['nit']])
            RV['funcalls'] = sp.array([info['funcalls']])
            self.null = RV
            if cache:
                f = h5py.File(out_file,'w')
                smartDumpDictHdf5(RV,f)
                f.close()
        return RV

    def getNull(self):
        """ get null model info """
        return self.null

    def setNull(self,null):
        """ set null model info """
        self.null = null

    ###########################################
    # Fitting alternative model
    ###########################################

    def optimize(self, G, params0=None, n_times=10, verbose=False, vmax=5, perturb=1e-3, factr=1e7):
        """
        Optimize the model considering G
        """
        # set params0 from null if params0 is None
        if params0 is None:
            if self.null is None:
                if verbose:     print(".. fitting null model")
                self.fitNull()
            if self.bgRE:
                params0 = sp.concatenate([self.null['params0_g'], self.null['params0_n']])
            else:
                params0 = self.null['params0_n']
            params_was_None = True
        else:
            params_was_None = False
        G *= sp.sqrt(self.N/(G**2).sum())
        self._gp.covar.G = G
        start = TIME.time()
        for i in range(n_times):
            if params_was_None:
                n_params = self.Cr.getNumberParams()
                _params0 = {'covar': sp.concatenate([1e-3*sp.randn(n_params), params0])}
            else:
                _params0 = {'covar': params0}
            self._gp.setParams(_params0)
            conv, info = self._gp.optimize(factr=factr, verbose=verbose)
            conv *= self.Cr.K().diagonal().max()<vmax
            conv *= self.getLMLgrad() < 0.1
            if conv or not params_was_None: break
        self.infoOpt = info
        if not conv:
            warnings.warn("not converged")
        # return value
        RV = {}
        if self.P>1:
            RV['Cr']  = self.Cr.K()
            if self.bgRE: RV['Cg']  = self.Cg.K()
            RV['Cn']  = self.Cn.K()
        RV['time']  = sp.array([TIME.time()-start])
        RV['params0'] = _params0
        RV['nit'] = sp.array([info['nit']])
        RV['funcalls'] = sp.array([info['funcalls']])
        RV['var']    = self.getVariances()
        RV['conv']  = sp.array([conv])
        RV['NLLAlt']  = sp.array([self.getNLLAlt()])
        RV['LLR']    = sp.array([self.getLLR()])
        RV['LMLgrad'] = sp.array([self.getLMLgrad()])
        return RV

    def getInfoOpt(self):
        """ get information for the optimization """
        return self.infoOpt

    def getVariances(self):
        """
        get variances
        """
        var = []
        var.append(self.Cr.K().diagonal())
        if self.bgRE:
            var.append(self.Cg.K().diagonal())
        var.append(self.Cn.K().diagonal())
        var = sp.array(var)
        return var

    def getNLLAlt(self):
        """
        get negative log likelihood of the alternative
        """
        return self._gp.LML()

    def getLLR(self):
        """
        get log likelihood ratio
        """
        assert self.null is not None, 'null model needs to be fitted!'
        return self.null['NLL0'][0] - self.getNLLAlt()

    def getLMLgrad(self):
        """
        get norm LML gradient
        """
        return (self._gp.LML_grad()['covar']**2).mean()

    def fitNullTraitByTrait(self, verbose=False, cache=False, out_dir='./cache', fname=None, rewrite=False):
        """
        Fit null model trait by trait
        """
        read_from_file = False
        if cache:
            assert fname is not None, 'MultiTraitSetTest:: specify fname'
            if not os.path.exists(out_dir): os.makedirs(out_dir)
            out_file = os.path.join(out_dir,fname)
            read_from_file = os.path.exists(out_file) and not rewrite

        RV = {}
        if read_from_file:
            f = h5py.File(out_file,'r')
            for p in range(self.P):
                trait_id = self.traitID[p]
                g = f[trait_id]
                RV[trait_id] = {}
                for key in list(g.keys()):
                    RV[trait_id][key] = g[key][:]
            f.close()
            self.nullST=RV
        else:
            """ create stSet and fit null column by column returns all info """
            if self.stSet is None:
                y = sp.zeros((self.N,1))
                self.stSet = MTSet(Y=y, S_R=self.S_R, U_R=self.U_R, F=self.F)
            RV = {}
            for p in range(self.P):
                trait_id = self.traitID[p]
                self.stSet.Y = self.Y[:,p:p+1]
                RV[trait_id] = self.stSet.fitNull()
            self.nullST = RV
            if cache:
                f = h5py.File(out_file,'w')
                smartDumpDictHdf5(RV,f)
                f.close()
        return RV

    def optimizeTraitByTrait(self, G, verbose=False, n_times=10, factr=1e3):
        """ Optimize trait by trait """
        assert self.nullST is not None, 'fit null model beforehand'
        RV = {}
        self.infoOptST = {}
        for p in range(self.P):
            trait_id = self.traitID[p]
            self.stSet.Y = self.Y[:, p:p+1]
            self.stSet.setNull(self.nullST[trait_id])
            RV[trait_id] = self.stSet.optimize(G, n_times=n_times, factr=factr, verbose=verbose)
            self.infoOptST[trait_id] = self.stSet.getInfoOpt()
        return RV

    def getInfoOptST(self):
        """ get information for the optimization """
        return self.infoOptST

    def _initParams(self, init_method=None):
        """ this function initializes the paramenter and Ifilter """
        if self.bgRE:
            if init_method=='random':
                params0 = {'covar': sp.randn(self._gpNull.covar.getNumberParams())}
            else:
                if self.P==1:
                    params0 = {'covar':sp.sqrt(0.5) * sp.ones(2)}
                else:
                    cov = 0.5*sp.cov(self.Y.T) + 1e-4*sp.eye(self.P)
                    chol = la.cholesky(cov, lower=True)
                    params = chol[sp.tril_indices(self.P)]
                    params0 = {'covar': sp.concatenate([params, params])}
        else:
            if self.P==1:
                params_cn = sp.array([1.])
            else:
                cov = sp.cov(self.Y.T) + 1e-4*sp.eye(self.P)
                chol = la.cholesky(cov, lower=True)
                params_cn = chol[sp.tril_indices(self.P)]
            params0 = {'covar': params_cn}
        return params0

if __name__=='__main__':
    from limix.utils.preprocess import covar_rescale

    sp.random.seed(0)

    #generate phenotype
    n = 1000
    p = 4
    f = 10
    Y = sp.randn(n, p)
    X = sp.randn(n, f)
    G = sp.randn(n, f)
    R = sp.dot(X, X.T)
    R = covar_rescale(R)

    mts = mtset(Y, R=R)
    nullMTInfo = mts.fitNull(cache=False)
    mts.optimize(G)
