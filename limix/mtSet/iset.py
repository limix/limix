import sys
import limix
import scipy as sp
import numpy as np
import numpy.linalg as nla
import scipy.linalg as la
from limix.mtSet.iset_full import ISet_Full
from limix.mtSet.iset_strat import ISet_Strat
import pandas as pd
from limix.mtSet.core.iset_utils import calc_emp_pv_eff
import pdb

def fit_iSet(Y, U_R=None, S_R=None, covs=None, Xr=None, n_perms=0, Ie=None,
             strat=False, verbose=True):
    """
    Args:
        Y:          [N, P] phenotype matrix
        S_R:        N vector of eigenvalues of R
        U_R:        [N, N] eigenvector matrix of R
        covs:       [N, K] matrix for K covariates
        Xr:         [N, S] genotype data of the set component
        n_perms:    number of permutations to consider
        Ie:         N boolean context indicator
        strat:      if True, the implementation with stratified designs is considered
    """
    factr=1e7 # remove?
    if strat:
        assert Ie is not None, 'Ie must be specified for stratification analyses'
        assert Y.shape[1]==1, 'Y must be Nx1 for stratification analysis'
    else:
        assert covs==None, 'Covariates are not supported for analysis of fully observed phenotypes'

    if verbose:     print('fittng iSet')
    if strat:
        mtSetGxE = ISet_Strat(Y, Ie, Xr, covs=covs)
        RV = {}
        RV['null'] = mtSetGxE.fitNull()
        RV['rank2'] = mtSetGxE.fitFullRank()
        RV['rank1'] = mtSetGxE.fitLowRank()
        RV['block'] = mtSetGxE.fitBlock()
        RV['var'] = mtSetGxE.getVC()
    else:
        mtSetGxE = ISet_Full(Y=Y, S_R=S_R, U_R=U_R, Xr=Xr, factr=factr)
        RV = {}
        RV['null'] = mtSetGxE.fitNull()
        RV['rank2'] = mtSetGxE.fitFullRank()
        RV['rank1'] = mtSetGxE.fitLowRank()
        LLR = RV['rank1']['NLLAlt'] - RV['rank2']['NLLAlt']
        if LLR<-1e-6:
            RV['rank2'] = mtSetGxE.fitFullRank(init_method='lr')
        try:
            RV['block'] = mtSetGxE.fitBlock()
        except:
            try:
                RV['block'] = mtSetGxE.fitBlock(init_method='null')
            except:
                RV['block'] = mtSetGxE.fitBlock(init_method='null_no_opt')
        RV['var'] = mtSetGxE.getVC()

    if n_perms>0:
        RVperm = {}
        nulls = ['null', 'block', 'rank1']
        tests = ['mtSet', 'iSet', 'iSet-het']
        for test in tests:
            RVperm[test+' LLR0'] = sp.zeros(n_perms)
        for seed_i in range(n_perms):
            if verbose:     print(('permutation %d / %d' % (seed_i, n_perms)))
            for it, test in enumerate(tests):
                if test=='mtSet':
                    idxs = sp.random.permutation(Xr.shape[0])
                    _Xr = Xr[idxs, :]
                    df0 = fit_iSet(Y, U_R=U_R, S_R=S_R, covs=covs, Xr=_Xr, n_perms=0, Ie=Ie, strat=strat, verbose=False)
                else:
                    Y0 = mtSetGxE._sim_from(set_covar=nulls[it])
                    Y0 -= Y0.mean(0)
                    df0 = fit_iSet(Y0, U_R=U_R, S_R=S_R, covs=covs, Xr=Xr, n_perms=0, Ie=Ie, strat=strat, verbose=False)
                RVperm[test+' LLR0'][seed_i]  = df0[test+' LLR'][0]

    # output
    LLR_mtSet = RV['null']['NLL']-RV['rank2']['NLL']
    LLR_iSet = RV['block']['NLL']-RV['rank2']['NLL']
    LLR_iSet_het = RV['rank1']['NLL']-RV['rank2']['NLL']

    if strat:   var_keys = ['var_r_full', 'var_c', 'var_n']
    else:       var_keys = ['var_r_full', 'var_g', 'var_n']
    varT = sp.sum([RV['var'][key] for key in var_keys])
    var_pers = RV['var']['var_r_block'] / varT
    var_resc = (RV['var']['var_r_rank1'] - RV['var']['var_r_block']) / varT
    var_het = (RV['var']['var_r_full'] - RV['var']['var_r_rank1']) / varT
    conv = RV['null']['conv']
    conv*= RV['block']['conv']
    conv*= RV['rank1']['conv']
    conv*= RV['rank2']['conv']

    M = sp.array([LLR_mtSet, LLR_iSet, LLR_iSet_het, var_pers, var_resc, var_het, conv]).T
    columns = ['mtSet LLR', 'iSet LLR', 'iSet-het LLR',
                'Persistent Var', 'Rescaling-GxC Var', 'Heterogeneity-GxC var', 'Converged']
    df = pd.DataFrame(M, columns=columns)

    if n_perms>0:
        return df, pd.DataFrame(RVperm)
    return df

if __name__=='__main__':

    N = 200
    C = 2
    S = 5
    K = 30
    Y = sp.randn(N,C)
    W = sp.randn(N,K)
    W-= W.mean(0)
    W/= W.std(0)
    R = sp.dot(W, W.T)
    R/= R.diagonal().mean(0)
    R+= 1e-4*sp.eye(N)
    S_R, U_R = nla.eigh(R)
    covs = U_R[:,-10:]
    Ie = sp.rand(200)<0.5

    df = pd.DataFrame()
    df0 = pd.DataFrame()
    n_regions = 5
    for i in range(n_regions):
        print(('.. analyzing region %d' % i))

        X = 1.*(sp.rand(N,S)<0.2)
        X-= X.mean(0)
        X/= X.std(0)
        X/= sp.sqrt(X.shape[1])
        #_df, _df0 = fit_iSet(Y, U_R, S_R, X, n_perms=10)
        _df, _df0 = fit_iSet(Y[:,[0]], Xr=X, covs=covs, n_perms=10, Ie=Ie, strat=True)
        df  = df.append(_df)
        df0 = df0.append(_df0)

    for test in ['mtSet', 'iSet', 'iSet-het']:
        df[test+' pv'] = calc_emp_pv_eff(df[test+' LLR'].values, df0[test+' LLR0'].values)
