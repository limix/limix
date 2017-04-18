import scipy as sp
import scipy.linalg as la
import pdb

def decompose_GxE(Cr):
    _S, _U = la.eigh(Cr)
    _w = sp.sqrt(_S[-1]) * _U[:,[-1]]
    RV = {}
    RV['r_full'] = Cr 
    RV['r_block'] = sp.mean(Cr) * sp.ones(_U.shape)
    RV['r_rank1'] = sp.dot(_w, _w.T)
    return RV

def msqrt(M):
    S, U = la.eigh(M)
    S[S<0] = 0
    Mh = U * S**(0.5)
    return Mh

def var_CoXX(C, X):
    C_h = msqrt(C)
    PK = sp.kron(C_h, X)
    PK-= PK.mean(0)
    rv = (PK**2).sum() / float(PK.shape[0] - 1.)
    return rv

def calc_emp_pv_eff(xt, x0):
    idxt = sp.argsort(xt)[::-1]
    idx0 = sp.argsort(x0)[::-1]
    xts = xt[idxt]
    x0s = x0[idx0]
    it = 0
    i0 = 0
    _count = 0
    count = sp.zeros(xt.shape[0])
    while 1:
        if x0s[i0]>xts[it]:
            _count += 1
            i0+=1
            if i0==x0.shape[0]:
                count[idxt[it:]] = _count
                break
        else:
            count[idxt[it]] = _count
            it+=1
            if it==xt.shape[0]:
                break
    RV = (count + 1) / float(x0.shape[0])
    RV[RV > 1.] = 1.
    return RV

