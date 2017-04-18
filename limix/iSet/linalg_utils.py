import scipy as sp
import scipy.linalg as la

def msqrt(C):
    _U, _S = la.eigh(C)
    _S[_S<0] = 0.
    return _U*sp.sqrt(_S)

def lowrank_approx(C, rank=1):
    _S, _U = la.eigh(C)
    _U = _U[:,-rank:]
    _S = _S[-rank:]
    R = sp.dot(_U*_S, _U.T)
    return R 

