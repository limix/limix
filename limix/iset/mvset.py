import limix
import scipy as sp
from .linalg_utils import lowrank_approx

ntype_dict = {'assoc': 'null', 'gxe': 'block', 'gxehet': 'rank1'}


def define_gp(Y, Xr, F, type, Rr):
    from limix_core.covar import LowRankCov
    from limix_core.covar import FixedCov
    from limix_core.covar import FreeFormCov
    from limix_core.gp import GP2KronSumLR
    from limix_core.gp import GP2KronSum

    P = Y.shape[1]
    _A = sp.eye(P)
    if type in ['null', 'rank1']:
        _Cr = LowRankCov(P, 1)
    elif type == 'block':
        _Cr = FixedCov(sp.ones((P, P)))
    elif type == 'full':
        _Cr = FreeFormCov(P)
    else:
        print('poppo')
    _Cn = FreeFormCov(P)
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
