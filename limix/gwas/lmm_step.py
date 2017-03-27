import scipy as sp
import scipy.stats as st
import scipy.linalg as la
import pdb
import limix
import time
from limix.core.gp import GP2KronSumLR
from limix.core.covar import FreeFormCov
from read_utils import read_geno

def calc_Ai_beta_s2(yKiy,FKiF,FKiy,df):
    Ai = la.pinv(FKiF)
    beta = sp.dot(Ai,FKiy)
    s2 = (yKiy-sp.dot(FKiy[:,0],beta[:,0]))/df
    return Ai,beta,s2

class LMMstep():

    def __init__(self,y,F,cov=None):
        if F is None:   F = sp.ones((y.shape[0],1))
        self.y = y
        self.F = F
        self.cov = cov
        self.df = y.shape[0]-F.shape[1]
        self._fit_null()

    def _fit_null(self):
        """ fit the null model """
        if self.cov==None:
            self.Kiy = self.y
            self.KiF = self.F
        else:
            self.Kiy = self.cov.solve(self.y)
            self.KiF = self.cov.solve(self.F)
        self.FKiy = sp.dot(self.F.T, self.Kiy)
        self.FKiF = sp.dot(self.F.T, self.KiF)
        self.yKiy = sp.dot(self.y[:,0], self.Kiy[:,0])
        # calc beta_F0 and s20
        self.A0i, self.beta_F0, self.s20 = calc_Ai_beta_s2(self.yKiy,self.FKiF,self.FKiy,self.df)

    def process(self, G, step=8, verbose=False):
        """ LMMstep scan """
        t0 = time.time()
        if (G.shape[1]%step!=0):
            'Error: step is not compatible with dim1(G)'
        n_tests = G.shape[1]/step
        # precompute some stuff
        if self.cov==None:  KiG = G
        else:               KiG = self.cov.solve(G)
        k = self.F.shape[1]
        GKiF = sp.dot(G.T,self.KiF)
        GKiy = sp.dot(G.T,self.Kiy)
        F1KiF1 = sp.zeros((k+step, k+step))
        F1KiF1[:k,:k] = self.FKiF
        F1Kiy = sp.zeros((k+step,1))
        F1Kiy[:k,0] = self.FKiy[:,0]
        s2 = sp.zeros(n_tests)
        self.beta_g = sp.zeros((step,n_tests))
        for i in range(n_tests):
            F1KiF1[k:,:k] = GKiF[i*step:(i+1)*step,:]
            F1KiF1[:k,k:] = F1KiF1[k:,:k].T
            F1KiF1[k:,k:] = sp.dot(G[:,i*step:(i+1)*step].T, KiG[:,i*step:(i+1)*step])
            F1Kiy[k:,0] = GKiy[i*step:(i+1)*step,0]
            _,beta,s2[i] = calc_Ai_beta_s2(self.yKiy,F1KiF1,F1Kiy,self.df)
            self.beta_g[:,i] = beta[k:,0]
        #dlml and pvs
        self.lrt = -self.df*sp.log(s2/self.s20)
        self.pv = st.chi2(1).sf(self.lrt)

        t1 = time.time()
        if verbose:
            print('Tested for %d variants in %.2f s' % (G.shape[1],t1-t0))

    def getPv(self):
        return self.pv

    def getBetaSNP(self):
        return self.beta_g

    def getBetaCov(self):
        return self.beta_F

    def getLRT(self):
        return self.lrt



if __name__=="__main__":

    N = 2000
    k = 10
    S = 1000
    y = sp.randn(N,1)
    E = sp.randn(N,k)
    G = 1.*(sp.rand(N,S)<0.2)
    F = sp.concatenate([sp.ones((N,1)), sp.randn(N,1)], 1)

    gp = GP2KronSumLR(Y=y, Cn=FreeFormCov(1), G=E, F=F, A=sp.ones((1,1)))
    gp.covar.Cr.setCovariance(0.5*sp.ones((1,1)))
    gp.covar.Cn.setCovariance(0.5*sp.ones((1,1)))
    gp.optimize()
    print('sg = %.2f' % gp.covar.Cr.K()[0,0])
    print('sn = %.2f' % gp.covar.Cn.K()[0,0])

    pdb.set_trace()

    print('New LMM')
    t0 = time.time()
    lmm = LMMstep(y,F,gp.covar)
    lmm.process(G)
    t1 = time.time()
    print('Elapsed:', t1-t0)
    pv = lmm.getPv()
    beta = lmm.getBetaSNP()
    lrt = lmm.getLRT()

    pdb.set_trace()
