# -*- coding: utf-8 -*-

import scipy as sp


class BF:
    r"""
    Calculates log Bayes factor between two models with different sets of environments included returning evidence for model 2 relative to model 1
    Parameters
    ----------
    y : (`N`, 1) ndarray
        phenotype vector
    x : (`N`, 1)
        SNP vector
    F : (`N`, L) ndarray
        fixed effect design for covariates.
    Env1 : (`N`, `K1`)
        Environmental matrix for model 1 (indviduals by number of environments to be included in model 1)
        If not set, the BF will be compare the environments included in model 2 to a model with no modelled GxE effects
    Env2 : (`N`, `K2`)
        Environmental matrix for model 2 (indviduals by number of environments to be included in model 2)
        If not set, the BF will be compare the environments included in model 1 to a model with no modelled GxE effects
    W : (`N`, `T`)
        design of random effect in the null model (identical for both models 1 and 2).
        By default, W is set to the larger of ``Env1``, ``Env2``.
    Examples
    --------
    This example shows how to run BF.
    .. doctest::
        >>> from numpy.random import RandomState
        >>> import scipy as sp
        >>> from limix.model.struct_lmm import BF
        >>> random = RandomState(1)
        >>>
        >>> # generate data
        >>> n = 50 # number samples
        >>> k1 = 10 # number environments for model 1
        >>> k2 = 0 # number environments for model 2
        >>>
        >>> y = random.randn(n, 1) # phenotype
        >>> x = 1. * (random.rand(n, 1) < 0.2) # genotype
        >>> E1 = random.randn(n, k1) # environemnts 1
        >>> E2 = random.randn(n, k2) # environemnts 1
        >>> covs = sp.ones((n, 1)) # intercept
        >>>
        >>> bf = BF(y, x, F = covs, Env1 = E1, Env2 = E2, W=E1)
        >>> bf.calc_bf()  # doctest: +FLOAT_CMP
        0.03013960889843048
    """

    def __init__(self, y, x, F, Env1, Env2, W=None):
        self.y = y
        self.x = x
        self.F = F
        self.Env1 = Env1
        self.Env2 = Env2
        self.W = W
        if self.W is None:
            if self.Env1.shape[1] <= self.Env2.shape[1]:
                self.W = self.Env2
            else:
                self.W = self.Env1

    def calc_lml(self, Env):
        from limix_core.covar import FreeFormCov
        from limix_core.gp import GP2KronSumLR

        _covs = sp.concatenate([self.F, self.W, self.x], 1)
        if Env.shape[1] == 0:
            xoE = sp.ones(self.x.shape)
        else:
            xoE = self.x * Env
        gp = GP2KronSumLR(Y=self.y, F=_covs, A=sp.eye(1), Cn=FreeFormCov(1), G=xoE)
        gp.covar.Cr.setCovariance(1e-4 * sp.ones((1, 1)))
        gp.covar.Cn.setCovariance(0.02 * sp.ones((1, 1)))
        gp.optimize(verbose=False)
        lml = -gp.LML()

        return lml

    def calc_bf(self):
        lml_model_1 = self.calc_lml(self.Env1)
        lml_model_2 = self.calc_lml(self.Env2)
        bf = lml_model_1 - lml_model_2
        return bf
