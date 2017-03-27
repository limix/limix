import sys
from limix.utils.preprocess import regressOut
import numpy as np

import scipy.linalg as LA
import copy

def compute_X1KX2(Y, D, X1, X2, A1=None, A2=None):

    R,C = Y.shape
    if A1 is None:
        nW_A1 = Y.shape[1]
        #A1 = np.eye(Y.shape[1])	#for now this creates A1 and A2
    else:
        nW_A1 = A1.shape[0]

    if A2 is None:
        nW_A2 = Y.shape[1]
        #A2 = np.eye(Y.shape[1])	#for now this creates A1 and A2
    else:
        nW_A2 = A2.shape[0]


    nW_X1 = X1.shape[1]
    rows_block = nW_A1 * nW_X1

    if 0:#independentX2:
        nW_X2 = 1
    else:
        nW_X2 = X2.shape[1]
    cols_block = nW_A2 * nW_X2

    block = np.zeros((rows_block,cols_block))


    if (R>C) or (A1 is None) or (A2 is None):
        for c in range(C):
            X1D = X1 * D[:,c:c+1]
            X1X2 = X1D.T.dot(X2)
            if (A1 is None) and (A2 is None):
                block[c*X1.shape[1]:(c+1)*X1.shape[1], c*X2.shape[1]:(c+1)*X2.shape[1]] += X1X2
            elif (A1 is None):
                block[c*X1.shape[1]:(c+1)*X1.shape[1],:] += np.kron(A2[:,c:c+1].T,X1X2)
            elif (A2 is None):
                block[:,c*X2.shape[1]:(c+1)*X2.shape[1]] += np.kron(A1[:,c:c+1],X1X2)
            else:
                A1A2 = np.outer(A1[:,c],A2[:,c])
                block += np.kron(A1A2,X1X2)
    else:
        for r in range(R):
            A1D = A1 * D[r:r+1,:]
            A1A2 = A1D.dot(A2.T)
            X1X2 = X1[r,:][:,np.newaxis].dot(X2[r,:][np.newaxis,:])
            block += np.kron(A1A2,X1X2)

    return block

class mean(cObject):

    def __init__(self,Y, identity_trick=False):
        """ init data term """
        self.Y = Y
        self.identity_trick=identity_trick
        self.clearFixedEffect()

    #########################################
    # Properties
    #########################################

    @property
    def A(self):
        return self._A

    @property
    def B(self):
        return self._B

    @property
    def F(self):
        return self._F

    @property
    def A_identity(self):
        return self._A_identity

    @property
    def REML_term(self):
        return self._REML_term

    @property
    def Y(self):
        return self._Y

    @property
    def N(self):
        return self._N

    @property
    def P(self):
        return self._P

    @property
    def n_fixed_effs(self):
        return self._n_fixed_effs

    @property
    def n_terms(self):
        return self._n_terms

    @property
    def Lr(self):
        return self._Lr

    @property
    def Lc(self):
        return self._Lc

    @property
    def d(self):
        return self._d

    @property
    def D(self):
        return np.reshape(self.d,(self.N,self.P), order='F')

    @property
    def LRLdiag(self):
        return self._LRLdiag

    @property
    def LCL(self):
        return self._LCL

    #########################################
    # Setters
    #########################################
    def use_identity_trick(self,identity_trick=True):
        self.identity_trick=identity_trick
        self.clear_cache('Fstar','Astar','Xstar','Xhat',
                         'Areml','Areml_eigh','Areml_chol','Areml_inv','beta_hat','B_hat',
                         'LRLdiag_Xhat_tens','Areml_grad',
                         'beta_grad','Xstar_beta_grad','Zstar','DLZ')


    def clearFixedEffect(self):
        """ erase all fixed effects """
        self._A = []
        self._F = []
        self._B = []
        self._A_identity = []
        self._REML_term = []
        self._n_terms = 0
        self._n_fixed_effs = 0
        self._n_fixed_effs_REML = 0
        self.indicator = {'term':np.array([]),
                            'row':np.array([]),
                            'col':np.array([])}
        self.clear_cache('Fstar','Astar','Xstar','Xhat',
                         'Areml','Areml_eigh','Areml_chol','Areml_inv','beta_hat','B_hat',
                         'LRLdiag_Xhat_tens','Areml_grad',
                         'beta_grad','Xstar_beta_grad','Zstar','DLZ')

    def addFixedEffect(self,F=None,A=None, REML=True, index=None):
        """
        set sample and trait designs
        F:      NxK sample design
        A:      LxP sample design
        REML:   REML for this term?
        index:  index of which fixed effect to replace. If None, just append.
        """
        if F is None:   F = np.ones((self.N,1))
        if A is None:
            A = np.eye(self.P)
            A_identity = True
        elif (A.shape == (self.P,self.P)) & (A==np.eye(self.P)).all():
            A_identity = True
        else:
            A_identity = False

        assert F.shape[0]==self.N, "F dimension mismatch"
        assert A.shape[1]==self.P, "A dimension mismatch"
        if index is None or index==self.n_terms:
            self.F.append(F)
            self.A.append(A)
            self.A_identity.append(A_identity)
            self.REML_term.append(REML)
            # build B matrix and indicator
            self.B.append(np.zeros((F.shape[1],A.shape[0])))
            self._n_terms+=1
            self._update_indicator(F.shape[1],A.shape[0])
        elif index >self.n_terms:
            raise Exception("index exceeds max index of terms")
        else:
            self._n_fixed_effs-=self.F[index].shape[1]*self.A[index].shape[0]
            if self.REML_term[index]:
                self._n_fixed_effs_REML-=self.F[index].shape[1]*self.A[index].shape[0]
            self.F[index] = F
            self.A[index] = A
            self.A_identity[index] = A_identity
            self.REML_term[index]=REML
            self.B[index] = np.zeros((F.shape[1],A.shape[0]))
            self._rebuild_indicator()

        self._n_fixed_effs+=F.shape[1]*A.shape[0]
        if REML:
            self._n_fixed_effs_REML+=F.shape[1]*A.shape[0]
        self.clear_cache('Fstar','Astar','Xstar','Xhat',
                         'Areml','Areml_eigh','Areml_chol','Areml_inv','beta_hat','B_hat',
                         'LRLdiag_Xhat_tens','Areml_grad',
                         'beta_grad','Xstar_beta_grad','Zstar','DLZ')

    def removeFixedEffect(self, index=None):
        """
        set sample and trait designs
        F:      NxK sample design
        A:      LxP sample design
        REML:   REML for this term?
        index:  index of which fixed effect to replace. If None, remove last term.
        """
        if self._n_terms==0:
            pass
        if index is None or index==(self._n_terms-1):

            self._n_terms-=1
            F = self._F.pop() #= self.F[:-1]
            A = self._A.pop() #= self.A[:-1]
            self._A_identity.pop() #= self.A_identity[:-1]
            REML_term = self._REML_term.pop()# = self.REML_term[:-1]
            self._B.pop()# = self.B[:-1]
            self._n_fixed_effs-=F.shape[1]*A.shape[0]
            if REML_term:
                self._n_fixed_effs_REML-=F.shape[1]*A.shape[0]

            pass
        elif index >= self.n_terms:
            raise Exception("index exceeds max index of terms")
        else:
            raise NotImplementedError("currently only last term can be removed")
            pass
        self._rebuild_indicator()
        self.clear_cache('Fstar','Astar','Xstar','Xhat',
                         'Areml','Areml_eigh','Areml_chol','Areml_inv','beta_hat','B_hat',
                         'LRLdiag_Xhat_tens','Areml_grad',
                         'beta_grad','Xstar_beta_grad','Zstar','DLZ')

    @Y.setter
    def Y(self,value):
        """ set phenotype """
        self._N,self._P = value.shape
        self._Y = value
        self.clear_cache('Ystar1','Ystar','Yhat','LRLdiag_Yhat',
                         'beta_grad','Xstar_beta_grad','Zstar','DLZ')

    @Lr.setter
    def Lr(self,value):
        """ set row rotation """
        assert value.shape[0]==self._N, 'dimension mismatch'
        assert value.shape[1]==self._N, 'dimension mismatch'
        self._Lr = value
        self.clear_cache('Fstar','Ystar1','Ystar','Yhat','Xstar','Xhat',
                         'Areml','Areml_eigh','Areml_chol','Areml_inv','beta_hat','B_hat',
                         'LRLdiag_Xhat_tens','LRLdiag_Yhat','Areml_grad',
                         'beta_grad','Xstar_beta_grad',
                         'LRLdiag_Xhat_tens','LRLdiag_Yhat','Areml_grad',
                         'beta_grad','Xstar_beta_grad','Zstar','DLZ')

    @Lc.setter
    def Lc(self,value):
        """ set col rotation """
        assert value.shape[0]==self._P, 'Lc dimension mismatch'
        assert value.shape[1]==self._P, 'Lc dimension mismatch'
        self._Lc = value
        self.clear_cache('Astar','Ystar','Yhat','Xstar','Xhat',
                         'Areml','Areml_eigh','Areml_chol','Areml_inv','beta_hat','B_hat',
                         'LRLdiag_Xhat_tens','LRLdiag_Yhat','Areml_grad',
                         'beta_grad','Xstar_beta_grad','Zstar','DLZ')

    @d.setter
    def d(self,value):
        """ set anisotropic scaling """
        assert value.shape[0]==self._P*self._N, 'd dimension mismatch'
        self._d = value
        self.clear_cache('Yhat','Xhat','Areml','Areml_eigh','Areml_chol','Areml_inv','beta_hat','B_hat',
                         'LRLdiag_Xhat_tens','LRLdiag_Yhat','Areml_grad',
                         'beta_grad','Xstar_beta_grad','Zstar','DLZ')

    @LRLdiag.setter
    def LRLdiag(self,value):
        """ set anisotropic scaling """
        self._LRLdiag = value
        self.clear_cache('LRLdiag_Xhat_tens','LRLdiag_Yhat','Areml_grad',
                         'beta_grad','Xstar_beta_grad')

    @LCL.setter
    def LCL(self,value):
        """ set anisotropic scaling """
        self._LCL = value
        self.clear_cache('Areml_grad','beta_grad','Xstar_beta_grad')

    #########################################
    # Getters (caching)
    #########################################

    @cached
    def Astar(self):
        RV = []
        for term_i in range(self.n_terms):
            RV.append(np.dot(self.A[term_i],self.Lc.T))
        return RV

    @cached
    def Fstar(self):
        RV = []
        for term_i in range(self.n_terms):
            RV.append(np.dot(self.Lr,self.F[term_i]))
        return RV

    def Ystar1(self):
        return np.dot(self.Lr,self.Y)

    @cached
    def Ystar(self):
        return np.dot(self.Ystar1(),self.Lc.T)

    @cached
    def Yhat(self):
        return self.D*self.Ystar()

    @cached
    def Xstar(self):
        RV = np.zeros((self.N*self.P,self.n_fixed_effs))
        ip = 0
        for i in range(self.n_terms):
            Ki = self.A[i].shape[0]*self.F[i].shape[1]
            RV[:,ip:ip+Ki] = np.kron(self.Astar()[i].T,self.Fstar()[i])
            ip += Ki
        return RV

    def var_total(self):
        return (self.Yhat()*self.Ystar()).sum()


    def var_explained(self):
        XKY = self.compute_XKY(M=self.Yhat())
        beta_hat = self.Areml_solve(XKY)
        return (XKY*beta_hat).sum(), beta_hat


    @cached
    def Xhat(self):
        RV = self.d[:,np.newaxis]*self.Xstar()
        return RV

    @cached
    def Areml(self):
        #A1 = self.XstarT_dot(self.Xhat())
        A2 =  self.compute_XKX()
        return A2

    @cached
    def Areml_chol(self):
        return LA.cholesky(self.Areml()).T

    @cached
    def Areml_REML_chol(self):
        return LA.cholesky(self.Areml()).T

    @cached
    def Areml_inv(self):
        return LA.cho_solve((self.Areml_chol(),True),np.eye(self.n_fixed_effs))

    #caching bug:
    #@cached
    def beta_hat(self):
        XKY = self.compute_XKY(M=self.Yhat())
        beta_hat = self.Areml_solve(XKY)
        return beta_hat


    @cached
    def B_hat(self):
        RV = []
        ip = 0
        for term_i in range(self.n_terms):
            RV.append(np.reshape(self.beta_hat()[ip:ip+self.B[term_i].size],self.B[term_i].shape, order='F'))
            ip += self.B[term_i].size
        return RV

    @cached
    def LRLdiag_Xhat_tens(self):
        RV  = np.reshape(self.Xhat(),(self.N,self.P,self.n_fixed_effs),order='F').copy()
        RV *= self.LRLdiag[:,np.newaxis,np.newaxis]
        return RV

    @cached
    def LRLdiag_Yhat(self):
        return self.LRLdiag[:,np.newaxis]*self.Yhat()

    @cached
    def Areml_grad(self):
        RV = np.einsum('jpk,lp->jlk',self.LRLdiag_Xhat_tens(),self.LCL)
        RV = RV.reshape((self.N*self.P,self.n_fixed_effs),order='F')
        RV*= self.d[:,np.newaxis]
        RV = -self.XstarT_dot(RV)
        return RV

    @cached
    def beta_grad(self):
        RV  = np.reshape(np.dot(self.LRLdiag_Yhat(),self.LCL.T),(self.N*self.P),order='F')
        RV *= self.d
        RV  = self.XstarT_dot(RV)
        RV += np.dot(self.Areml_grad(),self.beta_hat())
        RV  = -np.dot(self.Areml_inv(),RV)
        return RV

    @cached
    def Xstar_beta_grad(self):
        RV = np.zeros((self.N,self.P))
        ip = 0
        for term_i in range(self.n_terms):
            _Bgrad = np.reshape(self.beta_grad()[ip:ip+self.B[term_i].size],self.B[term_i].shape, order='F')
            RV+=np.dot(self.Fstar()[term_i],np.dot(_Bgrad,self.Astar()[term_i]))
            ip += self.B[term_i].size
        return RV


    @cached
    def Zstar(self):
        """ predict the value of the fixed effect """
        RV = self.Ystar().copy()
        for term_i in range(self.n_terms):
            if self.identity_trick and self.A_identity[term_i]:
                RV-=np.dot(self.Fstar()[term_i],self.B_hat()[term_i])
            else:
                RV-=np.dot(self.Fstar()[term_i],np.dot(self.B_hat()[term_i],self.Astar()[term_i]))
        self.clear_cache('DLZ')
        return RV

    @cached
    def Areml_eigh(self):
        """compute the eigenvalue decomposition of Astar"""
        s,U = LA.eigh(self.Areml(),lower=True)
        i_pos = (s>1e-10)
        s = s[i_pos]
        U = U[:,i_pos]
        return s,U

    @cached
    def DLZ(self):
        return self.Zstar()*np.reshape(self.D,(self.N,self.P), order='F')

    ###############################################
    # Other getters with no caching, should not they have caching somehow?
    ###############################################

    def Areml_solve(self, b):
        try:
            res = LA.cho_solve((self.Areml_chol(),True),b)
        except LA.LinAlgError:

            s,U = self.Areml_eigh()
            res = U.T.dot(b)
            res /= s[:,np.newaxis]
            res = U.dot(res)

        return res


    def compute_XKY(self, M=None):
        if M is None:
            M = self.Yhat()
        assert M.shape==(self.N,self.P)
        XKY = np.zeros((self.n_fixed_effs))
        n_weights = 0
        for term in range(self.n_terms):
            if self.identity_trick and self.A_identity[term]:
                XKY_block = compute_XYA(DY=M, X=self.Fstar()[term], A=None)
            else:
                XKY_block = compute_XYA(DY=M, X=self.Fstar()[term], A=self.Astar()[term])
            XKY[n_weights:n_weights + self.A[term].shape[0] * self.F[term].shape[1]] = XKY_block.ravel(order='F')
            n_weights += self.A[term].shape[0] * self.F[term].shape[1]
        return XKY

    def compute_XKX(self):
        #n_weights1 = 0
        #
        #for term1 in xrange(self.n_terms):
        #    n_weights1+=self.Astar()[term1].shape[0] * self.Fstar()[term1].shape[1]
        #cov_beta = np.zeros((n_weights1,n_weights1))
        cov_beta = np.zeros((self.n_fixed_effs,self.n_fixed_effs))
        n_weights1 = 0
        for term1 in range(self.n_terms):
            if self.identity_trick and self.A_identity[term1]:
                A_term1 = None
            else:
                A_term1 = self.Astar()[term1]
            n_weights2 = n_weights1
            for term2 in range(term1,self.n_terms):
                if self.identity_trick and self.A_identity[term2]:
                    A_term2 = None
                else:
                    A_term2 = self.Astar()[term2]
                block = compute_X1KX2(Y=self.Ystar(), D=self.D, X1=self.Fstar()[term1], X2=self.Fstar()[term2], A1=A_term1, A2=A_term2)
                cov_beta[n_weights1:n_weights1 + self.A[term1].shape[0] * self.F[term1].shape[1], n_weights2:n_weights2 + self.A[term2].shape[0] * self.F[term2].shape[1]] = block
                if term1!=term2:
                    cov_beta[n_weights2:n_weights2 + self.A[term2].shape[0] * self.F[term2].shape[1], n_weights1:n_weights1 + self.A[term1].shape[0] * self.F[term1].shape[1]] = block.T

                n_weights2+=self.A[term2].shape[0] * self.F[term2].shape[1]

            n_weights1+=self.A[term1].shape[0] * self.F[term1].shape[1]
        return cov_beta

    def predict(self):
        """ predict the value of the fixed effect """
        RV = np.zeros((self.N,self.P))
        for term_i in range(self.n_terms):
            RV+=np.dot(self.Fstar()[term_i],np.dot(self.B()[term_i],self.Astar()[term_i]))
        return RV

    def evaluate(self):
        """ predict the value of """
        RV  = -self.predict()
        RV += self.Ystar()
        return RV

    def getGradient(self,j):
        """ get rotated gradient for fixed effect i """
        i = int(self.indicator['term'][j])
        r = int(self.indicator['row'][j])
        c = int(self.indicator['col'][j])
        rv = -np.kron(self.Fstar()[i][:,[r]],self.Astar()[i][[c],:])
        return rv

    def XstarT_dot(self,M):
        """ get dot product of Xhat and M """
        if 0:
            #TODO: implement this properly
            pass
        else:
            RV = np.dot(self.Xstar().T,M)
        return RV

    def getResiduals(self):
        """ regress out fixed effects and results residuals """
        X = np.zeros((self.N*self.P,self.n_fixed_effs))
        ip = 0
        for i in range(self.n_terms):
            Ki = self.A[i].shape[0]*self.F[i].shape[1]
            X[:,ip:ip+Ki] = np.kron(self.A[i].T,self.F[i])
            ip += Ki
        y = np.reshape(self.Y,(self.Y.size,1),order='F')
        RV = regressOut(y,X)
        RV = np.reshape(RV,self.Y.shape,order='F')
        return RV

    #########################################
    # Params manipulation
    #########################################

    def getParams(self):
        """ get params """
        rv = np.array([])
        if self.n_terms>0:
            rv = np.concatenate([np.reshape(self.B[term_i],self.B[term_i].size, order='F') for term_i in range(self.n_terms)])
        return rv

    def setParams(self,params):
        """ set params """
        start = 0
        for i in range(self.n_terms):
            n_effects = self.B[i].size
            self.B[i] = np.reshape(params[start:start+n_effects],self.B[i].shape, order='F')
            start += n_effects

    #########################################
    # Utility functions
    #########################################

    def getDimensions(self):
        """ get phenotype dimensions """
        return self.N,self.P

    def _set_toChange(x):
        """ set variables in list x toChange """
        for key in list(x.keys()):
            self.toChange[key] = True

    def _update_indicator(self,K,L):
        """ update the indicator """
        _update = {'term': self.n_terms*np.ones((K,L)).T.ravel(),
                    'row': np.kron(np.arange(K)[:,np.newaxis],np.ones((1,L))).T.ravel(),
                    'col': np.kron(np.ones((K,1)),np.arange(L)[np.newaxis,:]).T.ravel()}
        for key in list(_update.keys()):
            self.indicator[key] = np.concatenate([self.indicator[key],_update[key]])

    def _rebuild_indicator(self):
        """ update the indicator """
        indicator = {'term':np.array([]),
                     'row':np.array([]),
                     'col':np.array([])}

        for term in range(self.n_terms):
            L = self.A[term].shape[0]
            K = self.F[term].shape[1]
            _update = {'term': (term+1)*np.ones((K,L)).T.ravel(),
                    'row': np.kron(np.arange(K)[:,np.newaxis],np.ones((1,L))).T.ravel(),
                    'col': np.kron(np.ones((K,1)),np.arange(L)[np.newaxis,:]).T.ravel()}
            for key in list(_update.keys()):
                indicator[key] = np.concatenate([indicator[key],_update[key]])
        self.indicator = indicator
