import sys
sys.path.append('./../../..')
import os
import subprocess
import pdb
import sys
import csv
import glob
import numpy as NP
from optparse import OptionParser
import time
import limix.stats.chi2mixture as C2M
import limix.utils.plot as plot
import scipy as SP

def plot_manhattan(pv,out_file):
    import matplotlib.pylab as PLT
    posCum = SP.arange(pv.shape[0])
    idx=~SP.isnan(pv[:,0])
    plot.plot_manhattan(posCum[idx],pv[idx][:,0],alphaNS=1.0,alphaS=1.0)
    PLT.savefig(out_file)


def postprocess(options):
    """ perform parametric fit of the test statistics and provide permutation and test pvalues """

    resdir = options.resdir
    out_file = options.outfile
    tol = options.tol

    print('.. load permutation results')
    file_name = os.path.join(resdir,'perm*','*.res')
    files = glob.glob(file_name)
    LLR0 = []
    for _file in files:
        print(_file)
        LLR0.append(NP.loadtxt(_file,usecols=[6]))
    LLR0 = NP.concatenate(LLR0)

    print('.. fit test statistics')
    t0 = time.time()
    c2m = C2M.Chi2mixture(tol=4e-3)
    c2m.estimate_chi2mixture(LLR0)
    pv0 = c2m.sf(LLR0)
    t1 = time.time()
    print(('finished in %s seconds'%(t1-t0)))

    print('.. export permutation results')
    perm_file = out_file+'.perm'
    RV = NP.array([LLR0,pv0]).T
    NP.savetxt(perm_file,RV,delimiter='\t',fmt='%.6f %.6e')

    print('.. load test results')
    file_name = os.path.join(resdir,'test','*.res')
    files = glob.glob(file_name)
    RV_test = []
    for _file in files:
        print(_file)
        RV_test.append(NP.loadtxt(_file))
    RV_test = NP.concatenate(RV_test)

    print('.. calc pvalues')
    pv = c2m.sf(RV_test[:,-1])[:,NP.newaxis]

    print('.. export test results')
    perm_file = out_file+'.test'
    RV_test = NP.hstack([RV_test,pv])
    NP.savetxt(perm_file,RV_test,delimiter='\t',fmt='%d %d %d %d %d %d %.6e %.6e')

    if options.manhattan:
        manhattan_file = out_file+'.manhattan.jpg'
        plot_manhattan(pv,manhattan_file)
