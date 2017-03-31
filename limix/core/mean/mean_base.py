import sys
from limix.core.type.observed import Observed
from limix.hcache import Cached, cached
from limix.utils.preprocess import regressOut
from limix.core.utils import assert_finite_array
import scipy as sp
import numpy as np

import scipy.linalg as LA
import copy

class MeanBase(Cached, Observed):
    """
    Basic mean function for gp regression
    Notation:
        N = number of individuals
        No = number of out-of-sample individuals for predictions
        K = number of fixed effect covariates
    """

    def __init__(self, y, W=None, Wstar=None):
        """
        Args:
            y:        phenotype vector [N, 1]
            W:        fixed effect design [N, K]
            Wstar:    out-of-sample fixed effect design for predictions [No, K]
        """
        Cached.__init__(self)
        self.y = y
        self.W = W
        self.Wstar = Wstar
        self.setFIinv(None)
        self._set_relay(None)

    def _set_relay(self, relay):
        self._relay = relay

    #########################################
    # Properties
    #########################################
    @property
    def W(self):
        return self._W

    @property
    def n_covs(self):
        return self.W.shape[1]

    @property
    def y(self):
        return self._y

    @property
    def b(self):
        return self._relay.b()

    @property
    def b_ste(self):
        if self.getFIinv() is None:
            R = None
        else:
            R = sp.sqrt(self.getFIinv().diagonal())[:,sp.newaxis]
        return R

    @property
    def n_data_points(self):
        return self._N

    # for multitrait interpretation
    @property
    def Y(self):
        return self.y

    @property
    def B(self):
        return self.b

    @property
    def B_ste(self):
        return self.b_ste

    @property
    def Wstar(self):
        return self._Wstar

    @property
    def use_to_predict(self):
        return self._use_to_predict

    #########################################
    # Setters
    #########################################
    #@Y.setter
    #def Y(self,value):
    #    """ set phenotype """
    #    self._N = value.shape[0]
    #    self._P = value.shape[1]
    #    self._Y = value
    #    self._notify()
    #    self.clear_cache('Yres')

    @W.setter
    def W(self,value):
        """ set fixed effect design """
        if value is None:   value = sp.zeros((self._N, 0))
        assert value.shape[0]==self._N, 'Dimension mismatch'
        self._K = value.shape[1]
        self._W = value
        self._notify()
        self.clear_cache('predict_in_sample','Yres')

    @Wstar.setter
    def Wstar(self,value):
        """ set fixed effect design for predictions """
        if value is None:
            self._use_to_predict = False
        else:
            assert value.shape[1]==self._K, 'Dimension mismatch'
            self._use_to_predict = True
        self._Wstar = value
        self.clear_cache('predict')

    @y.setter
    def y(self,value):
        assert_finite_array(value)
        assert value.shape[1] == 1, 'MeanBase: phenotype has to be a one column vector'
        self._N = value.shape[0]
        self._P = 1
        self._y = value
        # notify
        self._notify()
        self.clear_cache('Yres')

    @use_to_predict.setter
    def use_to_predict(self, value):
        assert not (self.Wstar is None and value is True), 'set Wstar!'
        self._use_to_predict = value

    #########################################
    # Predictions
    #########################################
    @cached
    def predict(self):
        r = self._predict_fun(self.Wstar)
        return r

    @cached
    def predict_in_sample(self):
        r = self._predict_fun(self.W)
        return r

    def _predict_fun(self,M):
        return sp.dot(M,self.b)

    def Yres(self):
        return self.yres

    @cached('Yres')
    def yres(self):
        return self.y - self.predict_in_sample()

    #######################################
    # Standard errors
    ########################################
    def setFIinv(self, value):
        self._FIinv = value

    def getFIinv(self):
        return self._FIinv

    ###########################################
    # Gradient TODO
    ###########################################
    #def getGradient(self,j):
    #    """ get gradient for fixed effect i """
    #    return rv

    #########################################
    # Params manipulation TODO
    #########################################
    #def getParams(self):
    #    """ get params """
    #    return rv

    #def setParams(self,params):
    #    """ set params """
    #    start = 0
    #    for i in range(self.n_terms):
    #        n_effects = self.B[i].size
    #        self.B[i] = np.reshape(params[start:start+n_effects],self.B[i].shape, order='F')
    #        start += n_effects
