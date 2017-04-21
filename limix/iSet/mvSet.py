import sys
import limix
from limix_core.covar import LowRankCov
from limix_core.covar import FixedCov
from limix_core.covar import FreeFormCov
from limix_core.gp import GP2KronSumLR
from limix_core.gp import GP2KronSum
import scipy as sp
import scipy.stats as st
from limix.mtSet.core.iset_utils import *
import numpy as np
import numpy.linalg as nla
import scipy.linalg as la
import copy
import pdb
from limix.util.preprocess import gaussianize
from scipy.optimize import fmin
import time
import pandas as pd
from .linalg_utils import msqrt
from .linalg_utils import lowrank_approx

ntype_dict = {'assoc': 'null', 'gxe': 'block', 'gxehet': 'rank1'}


def define_gp(Y, Xr, F, type, Rr):
    P = Y.shape[1]
    _A = sp.eye(P)
    if type in ['null', 'rank1']:
        _Cr = limix.core.covar.LowRankCov(P, 1)
    elif type == 'block':
        _Cr = limix.core.covar.FixedCov(sp.ones((P, P)))
    elif type == 'full':
        _Cr = limix.core.covar.FreeFormCov(P)
    else:
        print('poppo')
    _Cn = limix.core.covar.FreeFormCov(P)
    if type == 'null':
        _gp = GP2KronSumLR(
            Y=Y, G=sp.ones(
                (Y.shape[0], 1)), F=F, A=_A, Cr=_Cr, Cn=_Cn)
        _Cr.setParams(1e-9 * sp.ones(P))
        _gp.covar.act_Cr = False
    else:
        if Xr.shape[1] < Xr.shape[0]:
            _gp = GP2KronSumLR(Y=Y, G=Xr, F=F, A=_A, Cr=_Cr, Cn=_Cn)
        else:
            _gp = GP2KronSum(Y=Y, F=F, A=_A, R=Rr, Cg=_Cr, Cn=_Cn)
    return _gp


class MvSetTest():
    """
    Args:
        Y (ndarray):
            (`N`, `P`) phenotype matrix for `N` samples and `P` traits
        Xr (ndarray):
            (`N`, `S`) genotype values for `N` samples and `S` variants
            (defines the set component)
        factr (float):
            optimization paramenter that determines the accuracy of the solution
            (see scipy.optimize.fmin_l_bfgs_b for more details).
    """

    def __init__(
            self,
            Y=None,
            Xr=None,
            Rr=None,
            F=None,
            factr=1e7,
            debug=False):
        # avoid SVD failure by adding some jitter
        Xr += 2e-6 * (sp.rand(*Xr.shape) - 0.5)
        # make sure it is normalised
        Xr -= Xr.mean(0)
        Xr /= Xr.std(0)
        Xr /= sp.sqrt(Xr.shape[1])
        self.Y = Y
        self.F = F
        self.Xr = Xr
        self.covY = sp.cov(Y.T)
        self.factr = factr
        self.debug = debug
        self.gp = {}
        self.info = {}
        self.lowrank = Xr.shape[1] < Xr.shape[0]
        if Rr is not None:
            self.Rr = Rr
        if self.lowrank:
            self.Rr = None
        else:
            self.Rr = sp.dot(Xr, Xr.T)

    def assoc(self):
        # fit model
        for key in ['null', 'full']:
            if key not in list(self.gp.keys()):
                if self.debug:
                    print('.. dening %s' % key)
                self.gp[key] = define_gp(
                    self.Y, self.Xr, self.F, key, Rr=self.Rr)
                if self.debug:
                    print('.. fitting %s' % key)
                self.info[key] = self._fit(key, vc=True)
                # if key=='full':
                #    print self.info[key]['var'][0]
        return self.info['null']['LML'] - self.info['full']['LML']

    def gxe(self):
        # fit model
        for key in ['null', 'full', 'block']:
            if key not in list(self.gp.keys()):
                if self.debug:
                    print('.. defining %s' % key)
                self.gp[key] = define_gp(
                    self.Y, self.Xr, self.F, key, Rr=self.Rr)
                if self.debug:
                    print('.. fitting %s' % key)
                self.info[key] = self._fit(key, vc=True)
        return self.info['block']['LML'] - self.info['full']['LML']

    def gxehet(self):
        # fit model
        for key in ['null', 'full', 'rank1']:
            if key not in list(self.gp.keys()):
                if self.debug:
                    print('.. defining %s' % key)
                self.gp[key] = define_gp(
                    self.Y, self.Xr, self.F, key, Rr=self.Rr)
                if self.debug:
                    print('.. fitting %s' % key)
                self.info[key] = self._fit(key, vc=True)
        return self.info['rank1']['LML'] - self.info['full']['LML']

    def assoc_null(self, n_nulls=30):
        LLR0 = sp.zeros(n_nulls)
        for ni in range(n_nulls):
            idx_perms = sp.random.permutation(self.Y.shape[0])
            _Xr = self.Xr[idx_perms]
            if self.Rr is not None:
                _Rr = sp.dot(_Xr, _Xr.T)
            else:
                _Rr = None
            mvset0 = MvSetTest(Y=self.Y, F=self.F, Xr=_Xr, Rr=_Rr)
            LLR0[ni] = mvset0.assoc()
        return LLR0

    def gxe_null(self, n_nulls=30):
        LLR0 = sp.zeros(n_nulls)
        for ni in range(n_nulls):
            if self.lowrank:
                _Y = self.gp['block'].simulate_pheno()
            else:
                _Y = self.gp['block'].simulate_pheno(Rh=self.Xr)
            mvset0 = MvSetTest(Y=_Y, F=self.F, Xr=self.Xr, Rr=self.Rr)
            LLR0[ni] = mvset0.gxe()
        return LLR0

    def gxehet_null(self, n_nulls=30):
        LLR0 = sp.zeros(n_nulls)
        for ni in range(n_nulls):
            if self.lowrank:
                _Y = self.gp['rank1'].simulate_pheno()
            else:
                _Y = self.gp['rank1'].simulate_pheno(Rh=self.Xr)
            mvset0 = MvSetTest(Y=_Y, F=self.F, Xr=self.Xr, Rr=self.Rr)
            LLR0[ni] = mvset0.gxehet()
        return LLR0

    def _fit(self, type, vc=False):
        # 2. init
        if type == 'null':
            self.gp[type].covar.Cn.setCovariance(self.covY)
        elif type == 'full':
            Cn0_K = self.gp['null'].covar.Cn.K()
            # self.gp[type].covar.Cr.setCovariance(1e-4*sp.ones(self.covY.shape)+1e-4*sp.eye(self.covY.shape[0]))
            self.gp[type].covar.Cr.setCovariance(0.5 * Cn0_K)
            self.gp[type].covar.Cn.setCovariance(0.5 * Cn0_K)
        elif type == 'block':
            Crf_K = self.gp['full'].covar.Cr.K()
            Cnf_K = self.gp['full'].covar.Cn.K()
            self.gp[type].covar.Cr.scale = sp.mean(Crf_K)
            self.gp[type].covar.Cn.setCovariance(Cnf_K)
        elif type == 'rank1':
            Crf_K = self.gp['full'].covar.Cr.K()
            Cnf_K = self.gp['full'].covar.Cn.K()
            self.gp[type].covar.Cr.setCovariance(Crf_K)
            self.gp[type].covar.Cn.setCovariance(Cnf_K)
        else:
            print('poppo')
        self.gp[type].optimize(factr=self.factr, verbose=False)
        RV = {'Cr': self.gp[type].covar.Cr.K(),
              'Cn': self.gp[type].covar.Cn.K(),
              'B': self.gp[type].mean.B[0],
              'LML': sp.array([self.gp[type].LML()]),
              'LMLgrad': sp.array([sp.mean((self.gp[type].LML_grad()['covar'])**2)])}
        if vc:
            # tr(P CoR) = tr(C)tr(R) - tr(Ones C) tr(Ones R) / float(NP)
            #           = tr(C)tr(R) - C.sum() * R.sum() / float(NP)
            trRr = (self.Xr**2).sum()
            var_r = sp.trace(RV['Cr']) * trRr / float(self.Y.size - 1)
            var_c = sp.var(sp.dot(self.F, RV['B']))
            var_n = sp.trace(RV['Cn']) * self.Y.shape[0]
            var_n -= RV['Cn'].sum() / float(RV['Cn'].shape[0])
            var_n /= float(self.Y.size - 1)
            RV['var'] = sp.array([var_r, var_c, var_n])
            if 0 and self.Y.size < 5000:
                pdb.set_trace()
                Kr = sp.kron(RV['Cr'], sp.dot(self.Xr, self.Xr.T))
                Kn = sp.kron(RV['Cn'], sp.eye(self.Y.shape[0]))
                _var_r = sp.trace(Kr - Kr.mean(0)) / float(self.Y.size - 1)
                _var_n = sp.trace(Kn - Kn.mean(0)) / float(self.Y.size - 1)
                _var = sp.array([_var_r, var_c, _var_n])
                print(((_var - RV['var'])**2).mean())
            if type == 'full':
                # calculate within region vcs
                Cr_block = sp.mean(RV['Cr']) * sp.ones(RV['Cr'].shape)
                Cr_rank1 = lowrank_approx(RV['Cr'], rank=1)
                var_block = sp.trace(Cr_block) * trRr / float(self.Y.size - 1)
                var_rank1 = sp.trace(Cr_rank1) * trRr / float(self.Y.size - 1)
                RV['var_r'] = sp.array(
                    [var_block, var_rank1 - var_block, var_r - var_rank1])
        return RV


if 0:
    def _sim_from(self, set_covar='block', seed=None, qq=False):
        # 1. region term
        if set_covar == 'block':
            Cr = self.block['Cr']
            Cg = self.block['Cg']
            Cn = self.block['Cn']
        if set_covar == 'rank1':
            Cr = self.lr['Cr']
            Cg = self.lr['Cg']
            Cn = self.lr['Cn']
        Lc = msqrt(Cr)
        U, Sh, V = nla.svd(self.Xr, full_matrices=0)
        Lr = sp.zeros((self.Y.shape[0], self.Y.shape[0]))
        Lr[:, :Sh.shape[0]] = U * Sh[sp.newaxis, :]
        Z = sp.randn(*self.Y.shape)
        Yr = sp.dot(Lr, sp.dot(Z, Lc.T))
        # 2. bg term
        Lc = msqrt(Cg)
        Lr = self.XXh
        Z = sp.randn(*self.Y.shape)
        Yg = sp.dot(Lr, sp.dot(Z, Lc.T))
        # noise terms
        Lc = msqrt(Cn)
        Z = sp.randn(*self.Y.shape)
        Yn = sp.dot(Z, Lc.T)
        # normalize
        Y = Yr + Yg + Yn
        if qq:
            Y = gaussianize(Y)
            Y -= Y.mean(0)
            Y /= Y.std(0)
        return Y

if __name__ == '__main__':

    if 0:

        N = 10000
        P = 2
        S = 20
        Xr = 1. * (sp.rand(N, S) < 0.2)
        Y = sp.randn(N, P)
        F = sp.ones((N, 1))

        t0 = time.time()
        mvset = MvSetTest(Y=Y, Xr=Xr, F=F, factr=1e7)
        mvset.assoc()
        mvset.gxehet()
        print('.. permutations')
        mvset.assoc_null()
        print('.. bootstrap gxe')
        mvset.gxe_null()
        print('.. bootstrap gxehet')
        mvset.gxehet_null()
        print(time.time() - t0)

        pdb.set_trace()

    if 0:
        n_times = 100

        LLR_assoc = sp.zeros(n_times)
        LLR_gxe = sp.zeros(n_times)
        LLR_gxehet = sp.zeros(n_times)
        score_assoc = sp.zeros(n_times)
        score_gxe = sp.zeros(n_times)
        score_gxehet = sp.zeros(n_times)

        for time_i in range(n_times):
            print(time_i)

            N = 10000
            P = 2
            S = 20
            Xr = 1. * (sp.rand(N, S) < 0.2)
            Xr -= Xr.mean(0)
            Xr /= Xr.std(0)
            Xr /= sp.sqrt(Xr.shape[1])
            Y = sp.randn(N, P)
            F = sp.ones((N, 1))

            t0 = time.time()
            mvset = MvSetTest(Y=Y, Xr=Xr, F=F, factr=1e7)
            LLR_assoc[time_i] = mvset.assoc()
            LLR_gxe[time_i] = mvset.gxe()
            LLR_gxehet[time_i] = mvset.gxehet()
            print(time.time() - t0)

            t0 = time.time()
            score_assoc[time_i], _ = mvset.score(test='assoc')
            score_gxe[time_i], _ = mvset.score(test='gxe')
            score_gxehet[time_i], _ = mvset.score(test='gxehet')
            print(time.time() - t0)

        t0 = time.time()
        pdb.set_trace()
        import pylab as pl
        pl.ion()
        pl.subplot(221)
        pl.plot(LLR_assoc, score_assoc, '.k')
        pl.subplot(222)
        pl.plot(LLR_gxe, score_gxe, '.k')
        pl.subplot(223)
        pl.plot(LLR_gxehet, score_gxehet, '.k')
        pdb.set_trace()

    if 1:
        """
        Basic experiment for score test tests if optimisation of hyper matters
        """
        n_times = 100
        score_assoc = sp.zeros(n_times)
        score0_assoc = sp.zeros(n_times)
        score_gxehet = sp.zeros(n_times)
        score0_gxehet = sp.zeros(n_times)
        for time_i in range(n_times):
            print(time_i)
            N = 10000
            P = 2
            S = 20
            Xr = 1. * (sp.rand(N, S) < 0.2)
            Xr -= Xr.mean(0)
            Xr /= Xr.std(0)
            Xr /= sp.sqrt(Xr.shape[1])
            Y = sp.randn(N, P)
            F = sp.ones((N, 1))

            t0 = time.time()
            mvset = MvSetTest(Y=Y, Xr=Xr, F=F, factr=1e7)
            mvset.gxehet()
            score_assoc[time_i], _ = mvset.score(test='assoc')
            score0_assoc[time_i], _ = mvset.score(test='assoc', null=True)
            score_gxehet[time_i], _ = mvset.score(test='gxehet')
            score0_gxehet[time_i], _ = mvset.score(test='gxehet', null=True)

        pdb.set_trace()

    """
    Some stuff dor the Hessian that I am not using at the moment
    """
    if 0:
        # H = 0.5 * trace_term - quadratic_term
        # expected_H = -0.5 * trace_term
        # average_H = quadratic_term
        # 1. trace_term = tr(Ki_dKp_Ki_dKp)
        rankCr = 2
        T = []
        KiT = []
        for i in range(3):
            _S, _U = la.eigh(_Cr.K_grad_i(i))
            Cr_h = _U * sp.sqrt(_S)
            T.append(
                sp.reshape(
                    sp.kron(
                        Cr_h, Xr), (N, P, rankCr * S), order='F'))
            KiT.append(gp.covar.solve_t(T[i]))
        Ht = sp.zeros((3, 3))
        for i in range(3):
            Ht[i, i] = (sp.einsum('qpn,qpm->nm', T[i], KiT[i])**2).sum()
            for j in range(0, i):
                Ht[i, j] = (sp.einsum('qpn,qpm->nm', T[i], KiT[j])**2).sum()
                Ht[j, i] = Ht[i, j]
        if Y.shape[0] <= 1000:
            Ki = la.inv(
                sp.kron(
                    mvset.gp[type].covar.Cn.K(),
                    sp.eye(
                        Y.shape[0])))
            XrXr = sp.dot(Xr, Xr.T)
            KidK = [sp.dot(Ki, sp.kron(_Cr.K_grad_i(i), XrXr))
                    for i in range(3)]
            Ht0 = sp.zeros((3, 3))
            for i in range(3):
                for j in range(3):
                    Ht0[i, j] = sp.trace(sp.dot(KidK[i], KidK[j]))
            pdb.set_trace()

    if 0:
        # 2. quadratic_term = y_Ki_dKp_Ki_dKq_Ki_y
        Z = Y - sp.dot(F, mvset.gp[type].mean.B[0])
        XrXrKiZ = sp.dot(Xr, sp.dot(Xr.T, gp.covar.solve_t(Z)))
        XrXrKiZC = [sp.dot(XrXrKiZ, _Cr.K_grad_i(i).T) for i in range(3)]
        KiXrXrKiZC = [gp.covar.solve_t(XrXrKiZC[i]) for i in range(3)]
        Hq = sp.zeros((3, 3))
        for i in range(3):
            Hq[i, i] = (XrXrKiZC[i] * KiXrXrKiZC[i]).sum()
            for j in range(i):
                Hq[i, j] = (XrXrKiZC[i] * KiXrXrKiZC[j]).sum()
                Hq[j, i] = Hq[i, j]
        if Y.shape[0] <= 1000:
            z = sp.reshape(Z, (Z.size, 1), order='F')
            dK = [sp.kron(_Cr.K_grad_i(i), XrXr) for i in range(3)]
            Kiz = sp.dot(Ki, z)
            Hq0 = sp.zeros((3, 3))
            for i in range(3):
                for j in range(3):
                    Hq0[i, j] = sp.dot(z.T, sp.dot(
                        KidK[i], sp.dot(KidK[j], Kiz)))
            pdb.set_trace()

    if 0:
        # compute score with hessian

        H = 0.5 * Ht  # -0.5*Ht+Hq

        Hi = la.inv(H)
        score1[time_i] = (SSS * sp.dot(Hi, SSS)).sum()
        C = FreeFormCov(2)

        def f1(x):
            C.setParams(x)
            b = C.K()[sp.tril_indices(2)]
            delta = (b - SSS)
            val = (delta * sp.dot(Hi, delta)).sum()
            db_dx0 = C.K_grad_i(0)[sp.tril_indices(2)]
            db_dx1 = C.K_grad_i(1)[sp.tril_indices(2)]
            db_dx2 = C.K_grad_i(2)[sp.tril_indices(2)]
            grad = 2 * sp.array([(delta * sp.dot(Hi, db_dx0)).sum(),
                                 (delta * sp.dot(Hi, db_dx1)).sum(),
                                 (delta * sp.dot(Hi, db_dx2)).sum()])
            return val, grad
        x_opt, dscore1, info = sp.optimize.fmin_l_bfgs_b(f1, sp.randn(3))
        conv[time_i] = (info['grad']**2).mean() < 1e-5
        score1c[time_i] = score1[time_i] - dscore1
        print(time.time() - t0)
