"""Variance Decomposition testing code"""
import unittest
import scipy as SP
import scipy.stats
import pdb
import os
import sys
import limix
import limix.deprecated.modules.varianceDecomposition as VAR
from limix.test import data


class VarianceDecompoitionKronecker_test(unittest.TestCase):
    """test class for VarianceDecomposition
         THis testing class assumes everything is Kronecker.
    """

    def genGeno(self):
        X  = (SP.rand(self.N,self.S)<0.2)*1.
        self.D['X'] = X


    def genPheno(self):
        dp = SP.ones(self.P); dp[1]=-1
        Y = SP.zeros((self.N,self.P))
        gamma0_g = SP.randn(self.S)
        gamma0_n = SP.randn(self.N)
        for p in range(self.P):
            gamma_g = SP.randn(self.S)
            gamma_n = SP.randn(self.N)
            beta_g  = dp[p]*gamma0_g+gamma_g
            beta_n  = gamma0_n+gamma_n
            y=SP.dot(self.D['X'],beta_g)
            y+=beta_n
            Y[:,p]=y
        Y=SP.stats.zscore(Y,0)
        self.D['Y']= Y


    def setUp(self):
        #check: do we have a csv File?
        self.dir_name = os.path.dirname(__file__)
        self.dataset = os.path.join(self.dir_name,'varDecomp')

        if (not os.path.exists(self.dataset)) or 'recalc' in sys.argv:
            if not os.path.exists(self.dataset):
                os.makedirs(self.dataset)
            SP.random.seed(1)
            self.N = 200
            self.S = 1000
            self.P = 2
            self.D = {}
            self.genGeno()
            self.genPheno()
            self.generate = True
        else:
            self.generate=False
            #self.D = data.load(os.path.join(self.dir_name,self.dataset))
            self.D = data.load(self.dataset)
            self.N = self.D['X'].shape[0]
            self.S = self.D['X'].shape[1]
            self.P = self.D['Y'].shape[1]

        self.Kg = SP.dot(self.D['X'],self.D['X'].T)
        self.Kg = self.Kg/self.Kg.diagonal().mean()

        self.vc = VAR.VarianceDecomposition(self.D['Y'])
        self.vc.addRandomEffect(self.Kg,jitter=0)
        self.vc.addRandomEffect(is_noise=True,jitter=0)
        self.vc.addFixedEffect()

    def test_fit(self):
        """ optimization test """
        self.vc.optimize(verbose=False)
        params = self.vc.getScales()
        if self.generate:
            self.D['params_true'] = params
            data.dump(self.D,self.dataset)
            self.generate=False
        params_true = self.D['params_true']
        RV = ((SP.absolute(params)-SP.absolute(params_true))**2).max()
        self.assertTrue(RV<1e-6)

    def test_fitFast(self):
        """ optimization test """
        self.vc.optimize(fast=True,verbose=False)
        params = self.vc.getScales()
        if self.generate:
            self.D['params_true'] = params
            data.dump(self.D,self.dataset)
            self.generate=False

        params_true = self.D['params_true']
        #make sign invariant
        RV = ((SP.absolute(params)-SP.absolute(params_true))**2).max()<1e-6
        self.assertTrue(RV)

class VarianceDecompositionSoftKronecker_test(unittest.TestCase):
    """test class for VarianceDecomposition
         THis testing class assumes tests the ability to
    """

    def genGeno(self):
        X  = (SP.rand(self.N,self.S)<0.2)*1.
        self.D['X'] = X

    def genPheno(self):
        dp = SP.ones(self.P); dp[1]=-1
        Y = SP.zeros((self.N,self.P))
        gamma0_g = SP.randn(self.S)
        gamma0_n = SP.randn(self.N)
        for p in range(self.P):
            gamma_g = SP.randn(self.S)
            gamma_n = SP.randn(self.N)
            beta_g  = dp[p]*gamma0_g+gamma_g
            beta_n  = gamma0_n+gamma_n
            y=SP.dot(self.D['X'],beta_g)
            y+=beta_n
            Y[:,p]=y
        Y=SP.stats.zscore(Y,0)
        self.D['Y']= Y

    def setUp(self):
        #check: do we have a csv File?
        self.dir_name = os.path.dirname(__file__)
        self.dataset = os.path.join(self.dir_name,'varDecomp')

        if (not os.path.exists(self.dataset)) or 'recalc' in sys.argv:
            if not os.path.exists(self.dataset):
                os.makedirs(self.dataset)
            SP.random.seed(1)
            self.N = 200
            self.S = 1000
            self.P = 2
            self.D = {}
            self.genGeno()
            self.genPheno()
            self.generate = True
        else:
            self.generate=False
            #self.D = data.load(os.path.join(self.dir_name,self.dataset))
            self.D = data.load(self.dataset)
            self.N = self.D['X'].shape[0]
            self.S = self.D['X'].shape[1]
            self.P = self.D['Y'].shape[1]

        self.Kg = SP.dot(self.D['X'],self.D['X'].T)
        self.Kg = self.Kg/self.Kg.diagonal().mean()

        #add missing values to Y
        if 1:
            self.D['Y'][0,0] = SP.nan
            self.D['Y'][1,1] = SP.nan
            self.D['Y'][10,0] = SP.nan
            self.D['Y'][100,1] = SP.nan
        self.vc = VAR.VarianceDecomposition(self.D['Y'])
        self.vc.addRandomEffect(self.Kg,jitter=0)
        self.vc.addRandomEffect(is_noise=True,jitter=0)
        self.vc.addFixedEffect()

    def test_fit(self):
        """ optimization test """
        self.vc.optimize(verbose=False)
        params = self.vc.getScales()
        if self.generate:
            self.D['params_true'] = params
            data.dump(self.D,self.dataset)
            self.generate=False
        params_true = self.D['params_true']
        RV = ((SP.absolute(params)-SP.absolute(params_true))**2).max()
        #permit more flexibility, as we set a few values to NAN
        self.assertTrue(RV<1e-4)


if __name__ == '__main__':
    unittest.main()
