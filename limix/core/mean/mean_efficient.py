import sys
from limix.utils.preprocess import regressOut
from limix.core.utils.cobj import cObject
from limix.core.mean.mean import compute_X1KX2
from limix.core.mean.linear import compute_XYA
import  limix.utils.psd_solve as psd_solve
import numpy as np

import scipy.linalg as la
import copy


class Mean(cObject):

    def __init__(self,Y):
        """ init data term """
        self.Y = Y
        self.clearFixedEffect()

    #########################################
    # Properties
    #########################################

    @property
    def len(self):
        return len(self.A)

    @property
    def shape(self):
        return self.Y.shape

    def getDimensions(self):
        """ get phenotype dimensions """
        return self.shape

    @property
    def N(self):
        return self.Y.shape[0]

    @property
    def P(self):
        return self.Y.shape[1]

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

    @property
    def dof(self, index=None):
        """The number of degrees of freedom"""
        if index is None:
            dof = 0
            for i in range(self.len):
                dof += self.A[i].shape[0] * self.F[i].shape[1]
            return dof
        else:
            return self.A[index].shape[0] * self.F[index].shape[1]

    @property
    def dof_any(self):
        """The number of degrees of freedom for any effects terms"""
        return self.P * self.F_any.shape[1]

    @d.setter
    def d(self,value):
        """ set anisotropic scaling """
        assert value.shape[0]==self._P*self._N, 'd dimension mismatch'
        self._d = value
        self.clear_cache()
    #########################################
    # Setters
    #########################################

    def clearFixedEffect(self):
        """ erase all fixed effects """
        self.A = []
        self.F = []
        self.F_any = np.zeros((self.N,0))
        self.clear_cache()

    def addFixedEffect(self,F=None,A=None,index=None):
        """
        set sample and trait designs
        F:      NxK sample design
        A:      LxP sample design
        fast_computations:   False deactivates the fast computations for any and common effects (for debugging)
        """
        if F is None:
            F = np.ones((self.N,1))
        else:
            assert F.shape[0]==self.N, "F dimension mismatch"

        if ((A is None) or ( (A.shape == (self.P,self.P)) and (A==np.eye(self.P)).all() )):
            #case any effect
            self.F_any = np.hstack((self.F_any,F))
        elif (index is not None) and  ((A==self.A[index]).all()):
            #case common effect
            self.F[index] = np.hstack((self.F_index,F))
        else:
            #case general A
            assert A.shape[1]==self.P, "A dimension mismatch"
            self.F.append(F)
            self.A.append(A)

        self.clear_cache()


    @Lr.setter
    def Lr(self,value):
        """ set row rotation """
        assert value.shape==(self.N, self.N), 'dimension mismatch'
        self._Lr = value
        self.clear_cache()

    @Lc.setter
    def Lc(self,value):
        """ set col rotation """
        assert value.shape==(self.P, self.P), 'dimension mismatch'
        self._Lc = value
        self.clear_cache()

    @d.setter
    def d(self,value):
        """ set anisotropic scaling """
        assert value.shape[0]==self.P*self.N, 'd dimension mismatch'
        self._d = value
        self.clear_cache()

    @LRLdiag.setter
    def LRLdiag(self,value):
        """ set LRLdiag """
        self._LRLdiag = value
        self.clear_cache()

    @LCL.setter
    def LCL(self,value):
        """ set LCL """
        self._LCL = value
        self.clear_cache()

    #########################################
    # Getters (caching)
    #########################################

    @property
    def Astar(self):
        RV = []
        for term_i in range(self.len):
            RV.append(np.dot(self.A[term_i],self.Lc.T))
        return RV

    @property
    def Fstar(self):
        RV = []
        for term_i in range(self.len):
            RV.append(np.dot(self.Lr,self.F[term_i]))
        return RV

    @property
    def Astar_any(self):
        return np.eye(self.P)

    @property
    def Fstar_any(self):
        return np.dot(self.Lr,self.F_any)

    def Ystar1(self):
        return np.dot(self.Lr,self.Y)

    def Ystar(self):
        return np.dot(self.Ystar1(),self.Lc.T)

    def Yhat(self):
        return self.D*self.Ystar()

    @property
    def Areml_solver(self):
        return psd_solve.PsdSolverKron(A = self.XKX(), A_any = self.XanyKXany(), AA_any = self.XanyKX())

    def beta_hat(self):
        """compute ML beta"""
        XKY = self.XKY()
        XanyKY = self.XanyKY()
        beta_hat, beta_hat_any = self.Areml_solver.solve(b_any=XanyKY,b=XKY,check_finite=True)
        return beta_hat, beta_hat_any

    def var_total(self):
        return (self.Yhat()*self.Ystar()).sum()


    def var_explained(self, identity_trick = False):
        XKY = self.XKY()
        XanyKY = self.XanyKY()
        beta_hat, beta_hat_any = self.Areml_solver.solve(b_any=XanyKY,b=XKY,check_finite=True)
        var_explained = (XKY*beta_hat).sum() + (XanyKY*beta_hat_any).sum()
        return var_explained, beta_hat, beta_hat_any

    ###############################################
    # Other getters with no caching, should not they have caching somehow?
    ###############################################


    def XanyKXany(self):
        """
        compute self covariance for any
        """
        result = np.empty((self.P,self.F_any.shape[1],self.F_any.shape[1]), order='C')
        for p in range(self.P):
            X1D = self.Fstar_any * self.D[:,p:p+1]
            X1X2 = X1D.T.dot(self.Fstar_any)
            result[p] = X1X2
        return result


    def XanyKX(self):
        """
        compute cross covariance for any and rest
        """
        result = np.empty((self.P,self.F_any.shape[1],self.dof), order='C')
        #This is trivially parallelizable:
        for p in range(self.P):
            FanyD = self.Fstar_any * self.D[:,p:p+1]
            start = 0
            #This is trivially parallelizable:
            for term in range(self.len):
                stop = start + self.F[term].shape[1]*self.A[term].shape[0]
                result[p,:,start:stop] = self.XanyKX2_single_p_single_term(p=p, F1=FanyD, F2=self.Fstar[term], A2=self.Astar[term])
                start = stop
        return result

    def XKX(self):
        """
        compute self covariance for rest
        """
        cov_beta = np.zeros((self.dof,self.dof))
        start_row = 0
        #This is trivially parallelizable:
        for term1 in range(self.len):
            stop_row = start_row + self.A[term1].shape[0] * self.F[term1].shape[1]
            start_col = start_row
            #This is trivially parallelizable:
            for term2 in range(term1,self.len):
                stop_col = start_col + self.A[term2].shape[0] * self.F[term2].shape[1]
                cov_beta[start_row:stop_row, start_col:stop_col] = compute_X1KX2(Y=self.Ystar(), D=self.D, X1=self.Fstar[term1], X2=self.Fstar[term2], A1=self.Astar[term1], A2=self.Astar[term2])
                if term1!=term2:
                    cov_beta[start_col:stop_col, start_row:stop_row] = cov_beta[n_weights1:stop_row, n_weights2:stop_col].T
                start_col = stop_col
            start_row = stop_row
        return cov_beta


    def XanyKX2_single_p_single_term(self, p, F2, A2=None, F1=None):
            X1D = F1
            if F1 is None:
                self.Fstar_any * self.D[:,p:p+1]
            X1X2 = X1D.T.dot(F2)
            if A2 is not None:
                return np.kron(A2[:,p:p+1].T,X1X2)
            else:
                return X1X2

    def XanyKY(self):

        return compute_XYA(DY=self.Yhat(), X=self.Fstar_any, A=None).T

    def XKY(self):
        M = self.Yhat()
        XKY = np.zeros((self.dof))

        n_weights = 0
        for term in range(self.len):
            XKY_block = compute_XYA(DY=M, X=self.Fstar[term], A=self.Astar[term])
            XKY[n_weights:n_weights + self.A[term].shape[0] * self.F[term].shape[1]] = XKY_block.ravel(order='F')
            n_weights += self.A[term].shape[0] * self.F[term].shape[1]

        return XKY


    #########################################
    # Utility functions
    #########################################
