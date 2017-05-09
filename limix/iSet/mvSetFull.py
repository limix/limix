import scipy as sp
import scipy.linalg as la
from .linalg_utils import lowrank_approx

ntype_dict = {'assoc': 'null', 'gxe': 'block', 'gxehet': 'rank1'}


def define_gp(Y, Xr, Sg, Ug, type):
    from limix_core.covar import LowRankCov
    from limix_core.covar import FixedCov
    from limix_core.covar import FreeFormCov
    from limix_core.gp import GP3KronSumLR
    from limix_core.gp import GP2KronSum
    P = Y.shape[1]
    _A = sp.eye(P)
    if type in 'rank1':
        _Cr = LowRankCov(P, 1)
    elif type == 'block':
        _Cr = FixedCov(sp.ones((P, P)))
    elif type == 'full':
        _Cr = FreeFormCov(P)
    elif type == 'null':
        pass
    else:
        print('poppo')
    _Cn = FreeFormCov(P)
    _Cg = FreeFormCov(P)
    if type == 'null':
        _gp = GP2KronSum(Y=Y, Cg=_Cg, Cn=_Cn, S_R=Sg, U_R=Ug)
    else:
        _gp = GP3KronSumLR(Y=Y, G=Xr, Cr=_Cr, Cg=_Cg, Cn=_Cn, S_R=Sg, U_R=Ug)
    return _gp


class MvSetTestFull():

    def __init__(
            self,
            Y=None,
            Xr=None,
            Rg=None,
            Ug=None,
            Sg=None,
            factr=1e7,
            debug=False):
        """
        Args:
            Y:          [N, P] phenotype matrix
            Xr:         [N, S] genotype data of the set component
            R:          [N, S] genotype data of the set component
            factr:      paramenter that determines the accuracy of the solution
                        (see scipy.optimize.fmin_l_bfgs_b for more details)
        """
        # assert Xr
        Xr -= Xr.mean(0)
        Xr /= Xr.std(0)
        Xr /= sp.sqrt(Xr.shape[1])
        self.Y = Y
        self.Xr = Xr
        if Sg is None or Ug is None:
            Sg, Ug = la.eigh(Rg)
        self.Rg = Rg
        self.Ug = Ug
        self.Sg = Sg
        self.covY = sp.cov(Y.T)
        self.factr = factr
        self.debug = debug
        self.gp = {}
        self.info = {}
        # _trRr = sp.diagonal(sp.dot(self.Ug, sp.dot(sp.diag(self.Sg), self.Ug.T))).sum()
        self.trRg = ((self.Ug * self.Sg**0.5)**2).sum()

    def assoc(self):
        # fit model
        for key in ['null', 'full']:
            if key not in list(self.gp.keys()):
                if self.debug:
                    print('.. dening %s' % key)
                self.gp[key] = define_gp(
                    self.Y, self.Xr, self.Sg, self.Ug, key)
                if self.debug:
                    print('.. fitting %s' % key)
                self.info[key] = self._fit(key, vc=True)
        return self.info['null']['LML'] - self.info['full']['LML']

    def gxe(self):
        # fit model
        for key in ['null', 'full', 'block']:
            if key not in list(self.gp.keys()):
                if self.debug:
                    print('.. defining %s' % key)
                self.gp[key] = define_gp(
                    self.Y, self.Xr, self.Sg, self.Ug, key)
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
                    self.Y, self.Xr, self.Sg, self.Ug, key)
                if self.debug:
                    print('.. fitting %s' % key)
                self.info[key] = self._fit(key, vc=True)
        return self.info['rank1']['LML'] - self.info['full']['LML']

    def assoc_null(self, n_nulls=30):
        LLR0 = sp.zeros(n_nulls)
        for ni in range(n_nulls):
            idx_perms = sp.random.permutation(self.Y.shape[0])
            _Xr = self.Xr[idx_perms]
            mvset0 = MvSetTestFull(Y=self.Y, Sg=self.Sg, Ug=self.Ug, Xr=_Xr)
            LLR0[ni] = mvset0.assoc()
        return LLR0

    def gxe_null(self, n_nulls=30):
        LLR0 = sp.zeros(n_nulls)
        for ni in range(n_nulls):
            _Y = self.gp['block'].simulate_pheno()
            mvset0 = MvSetTestFull(Y=_Y, Sg=self.Sg, Ug=self.Ug, Xr=self.Xr)
            LLR0[ni] = mvset0.gxe()
        return LLR0

    def gxehet_null(self, n_nulls=30):
        LLR0 = sp.zeros(n_nulls)
        for ni in range(n_nulls):
            _Y = self.gp['rank1'].simulate_pheno()
            mvset0 = MvSetTestFull(Y=_Y, Sg=self.Sg, Ug=self.Ug, Xr=self.Xr)
            LLR0[ni] = mvset0.gxehet()
        return LLR0

    def score(self):
        pass

    def _fit(self, type, vc=False):
        # 2. init
        if type == 'null':
            self.gp[type].covar.Cg.setCovariance(0.5 * self.covY)
            self.gp[type].covar.Cn.setCovariance(0.5 * self.covY)
        elif type == 'full':
            Cg0_K = self.gp['null'].covar.Cg.K()
            Cn0_K = self.gp['null'].covar.Cn.K()
            self.gp[type].covar.Cr.setCovariance((Cn0_K + Cg0_K) / 3.)
            self.gp[type].covar.Cg.setCovariance(2. * Cg0_K / 3.)
            self.gp[type].covar.Cn.setCovariance(2. * Cn0_K / 3.)
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
        RV = {'Cg': self.gp[type].covar.Cg.K(),
              'Cn': self.gp[type].covar.Cn.K(),
              'LML': sp.array([self.gp[type].LML()]),
              'LMLgrad': sp.array([sp.mean((self.gp[type].LML_grad()['covar'])**2)])}
        if type == 'null':
            RV['Cr'] = sp.zeros(RV['Cg'].shape)
        else:
            RV['Cr'] = self.gp[type].covar.Cr.K()
        if vc:
            # tr(P CoR) = tr(C)tr(R) - tr(Ones C) tr(Ones R) / float(NP)
            #           = tr(C)tr(R) - C.sum() * R.sum() / float(NP)
            trRr = (self.Xr**2).sum()
            var_r = sp.trace(RV['Cr']) * trRr / float(self.Y.size - 1)
            var_g = sp.trace(RV['Cg']) * self.trRg / float(self.Y.size - 1)
            var_n = sp.trace(RV['Cn']) * self.Y.shape[0]
            var_n -= RV['Cn'].sum() / float(RV['Cn'].shape[0])
            var_n /= float(self.Y.size - 1)
            RV['var'] = sp.array([var_r, var_g, var_n])
            if 0 and self.Y.size < 5000:
                pdb.set_trace()
                Kr = sp.kron(RV['Cr'], sp.dot(self.Xr, self.Xr.T))
                Kn = sp.kron(RV['Cn'], sp.eye(self.Y.shape[0]))
                _var_r = sp.trace(Kr - Kr.mean(0)) / float(self.Y.size - 1)
                _var_n = sp.trace(Kn - Kn.mean(0)) / float(self.Y.size - 1)
                _var = sp.array([_var_r, var_g, _var_n])
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
