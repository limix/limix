import scipy as sp
import scipy.stats as st
import scipy.linalg as la
import pdb
import limix
import time
from limix.core.gp import GP2KronSumLR
from limix.core.covar import FreeFormCov
from read_utils import read_geno

class LMM():

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
        self.A0i = la.inv(self.FKiF)
        self.beta_F0 = sp.dot(self.A0i,self.FKiy)
        self.s20 = (self.yKiy-sp.dot(self.FKiy[:,0],self.beta_F0[:,0]))/self.df 

    def process(self, G, verbose=False):
        """ LMM scan """
        t0 = time.time()
        # precompute some stuff
        if self.cov==None:  KiG = G
        else:               KiG = self.cov.solve(G)
        GKiy = sp.dot(G.T, self.Kiy[:,0])
        GKiG = sp.einsum('ij,ij->j', G, KiG)
        FKiG = sp.dot(self.F.T, KiG)

        # Let us denote the inverse of Areml as 
        # Ainv = [[A0i + m mt / n, m], [mT, n]]
        A0iFKiG = sp.dot(self.A0i, FKiG) 
        n = 1./(GKiG-sp.einsum('ij,ij->j', FKiG, A0iFKiG))
        M = -n*A0iFKiG
        self.beta_F = self.beta_F0+M*sp.dot(M.T,self.FKiy[:,0])/n
        self.beta_F+= M*GKiy 
        self.beta_g = sp.einsum('is,i->s',M,self.FKiy[:,0])
        self.beta_g+= n*GKiy 

        # sigma
        s2 = self.yKiy-sp.einsum('i,is->s',self.FKiy[:,0],self.beta_F)
        s2-= GKiy*self.beta_g
        s2/= self.df
        
        #dlml and pvs
        self.lrt = -self.df*sp.log(s2/self.s20)
        self.pv = st.chi2(1).sf(self.lrt)

        t1 = time.time()
        if verbose:
            print 'Tested for %d variants in %.2f s' % (G.shape[1],t1-t0)

    def process_in_batch(self, geno, block_size=5000, type='bed', verbose=False):
        t0 = time.time()
        idxs = sp.arange(0,geno.shape[0],block_size)
        n_blocks = idxs.shape[0]
        idxs = sp.append(idxs, geno.shape[0])
        beta_F = sp.zeros((self.F.shape[1], geno.shape[0]))
        beta_g = sp.zeros(geno.shape[0])
        lrt = sp.zeros(geno.shape[0]) 
        pv = sp.zeros(geno.shape[0]) 
        for block_i in range(n_blocks):
            print 'Block %d/%d' % (block_i, n_blocks)
            idx0 = idxs[block_i]
            idx1 = idxs[block_i+1]
            _G = read_geno(geno,idx0,idx1,type)
            self.process(_G)
            beta_F[:,idx0:idx1] = self.getBetaCov()
            beta_g[idx0:idx1] = self.getBetaSNP()
            lrt[idx0:idx1] = self.getLRT()
            pv[idx0:idx1] = self.getPv()
        self.beta_F = beta_F
        self.beta_g = beta_g
        self.lrt = lrt
        self.pv = pv
        t1 = time.time()
        if verbose:
            print 'Tested for %d variants in %.2f s' % (pv.shape[0],t1-t0)
            
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
    print 'sg = %.2f' % gp.covar.Cr.K()[0,0]
    print 'sn = %.2f' % gp.covar.Cn.K()[0,0]

    print 'New LMM'
    t0 = time.time()
    lmm = LMM(y,F,gp.covar)
    lmm.process(G)
    t1 = time.time()
    print 'Elapsed:', t1-t0
    pv = lmm.getPv()
    beta = lmm.getBetaSNP()
    lrt = lmm.getLRT()

