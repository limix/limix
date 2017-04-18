"""@package fdr
FDR estimation using Benjamini Hochberg and Stories method
"""

import numpy as np
import scipy as sp
import sys, pickle, pdb
import scipy.stats as st
import scipy.interpolate
import logging as LG


def qvalues1(PV,m=None,pi=1.0):
    """estimate q vlaues from a list of Pvalues
    this algorihm is taken from Storey, significance testing for genomic ...
    m: number of tests, (if not len(PV)), pi: fraction of expected true null (1.0 is a conservative estimate)
    @param PV: pvalues
    @param m:  total number of tests if PV is not the entire array.
    @param pi: fraction of true null
    """

    S = PV.shape
    PV = PV.flatten()
    if m is None:
        m = len(PV) * 1.0
    else:
        m*=1.0
    lPV = len(PV)

    #1. sort pvalues
    PV = PV.squeeze()
    IPV = PV.argsort()
    PV  = PV[IPV]

    #2. estimate lambda
    if pi is None:
        lrange = sp.linspace(0.05,0.95,max(lPV/100.0,10))
        pil    = sp.double((PV[:,sp.newaxis]>lrange).sum(axis=0))/lPV
        pilr   = pil/(1.0-lrange)
        #ok, I think for SNPs this is pretty useless, pi is close to 1!
        pi =1.0
        #if there is something useful in there use the something close to 1
        if pilr[-1]<1.0:
            pi = pilr[-1]

    #3. initialise q values
    QV_ = pi * m/lPV* PV
    QV_[-1] = min(QV_[-1],1.0)
    #4. update estimate
    for i in range(lPV-2,-1,-1):
        QV_[i] = min(pi*m*PV[i]/(i+1.0),QV_[i+1])
    #5. invert sorting
    QV = sp.zeros_like(PV)
    QV[IPV] = QV_

    QV = QV.reshape(S)
    return QV

def qvalues(pv, m = None, return_pi0 = False, lowmem = False, pi0 = None, fix_lambda = None):

    original_shape = pv.shape

    assert(pv.min() >= 0 and pv.max() <= 1)

    pv = pv.ravel() # flattens the array in place, more efficient than flatten()

    if m == None:
        m = float(len(pv))
    else:
        # the user has supplied an m, let's use it
        m *= 1.0

    # if the number of hypotheses is small, just set pi0 to 1
    if len(pv) < 100:
        pi0 = 1.0
    elif pi0 != None:
        pi0 = pi0
    else:
        # evaluate pi0 for different lambdas
        pi0 = []
        lam = sp.arange(0, 0.90, 0.01)
        counts = sp.array([(pv > i).sum() for i in sp.arange(0, 0.9, 0.01)])

        if fix_lambda != None:
            interv_count = (pv > fix_lambda - 0.01).sum()
            uniform_sim = sp.array([(pv > fix_lambda-0.01).sum()*(i+1) for i in sp.arange(0, len(sp.arange(0, 0.90, 0.01)))][::-1])
            counts += uniform_sim

        for l in range(len(lam)):
            pi0.append(counts[l]/(m*(1-lam[l])))

        pi0 = sp.array(pi0)

        # fit natural cubic spline
        tck = sp.interpolate.splrep(lam, pi0, k = 3)
        pi0 = sp.interpolate.splev(lam[-1], tck)
        if pi0 > 1:
            LG.warning("got pi0 > 1 (%.3f) while estimating qvalues, setting it to 1" % pi0)
            pi0 = 1.0

        assert(pi0 >= 0 and pi0 <= 1), "%f" % pi0


    if lowmem:
        # low memory version, only uses 1 pv and 1 qv matrices
        qv = sp.zeros((len(pv),))
        last_pv = pv.argmax()
        qv[last_pv] = (pi0*pv[last_pv]*m)/float(m)
        pv[last_pv] = -sp.inf
        prev_qv = last_pv
        for i in range(int(len(pv))-2, -1, -1):
            cur_max = pv.argmax()
            qv_i = (pi0*m*pv[cur_max]/float(i+1))
            pv[cur_max] = -sp.inf
            qv_i1 = prev_qv
            qv[cur_max] = min(qv_i, qv_i1)
            prev_qv = qv[cur_max]

    else:
        p_ordered = sp.argsort(pv)
        pv = pv[p_ordered]
        # estimate qvalues
    #     qv = pi0*m*pv/(sp.arange(len(pv))+1.0)

    #     for i in xrange(int(len(qv))-2, 0, -1):
    #         qv[i] = min([qv[i], qv[i+1]])


    qv = pi0 * m/len(pv) * pv
    qv[-1] = min(qv[-1],1.0)

    for i in range(len(pv)-2, -1, -1):
        qv[i] = min(pi0*m*pv[i]/(i+1.0), qv[i+1])



    # reorder qvalues
    qv_temp = qv.copy()
    qv = sp.zeros_like(qv)
    qv[p_ordered] = qv_temp



    # reshape qvalues
    qv = qv.reshape(original_shape)

    if return_pi0:
        return qv, pi0
    else:
        return qv


def estimate_lambda(pv):
    """estimate lambda form a set of PV"""
    LOD2 = sp.median(st.chi2.isf(pv,1))
    L = (LOD2/0.456)
    return (L)

def LOD2PV(lods):
    PV = (st.chi2.sf(2*lods, 1))
    return PV
