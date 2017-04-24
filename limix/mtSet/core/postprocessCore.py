import sys
import os
import subprocess
import pdb
import sys
import csv
import glob
import numpy as NP
from optparse import OptionParser
import time
import scipy as SP
import pandas as pd


def postprocess(options):
    r"""
    perform parametric fit of the test statistics and provide
    permutation and test pvalues
    """

    from limix.stats import Chi2mixture

    import pdb; pdb.set_trace()

    resdir = options.resdir
    out_file = options.outfile
    tol = options.tol

    print('.. load permutation results')
    file_name = os.path.join(resdir, 'perm*', '*.res')
    files = glob.glob(file_name)
    LLR0 = []
    for _file in files:
        print(_file)
        _LLR0 = pd.DataFrame.from_csv(_file, sep='\t')['LLR'].values
        LLR0.append(_LLR0)
    LLR0 = NP.concatenate(LLR0)

    print('.. fit test statistics')
    t0 = time.time()
    c2m = Chi2mixture(tol=4e-3)
    c2m.estimate_chi2mixture(LLR0)
    pv0 = c2m.sf(LLR0)
    t1 = time.time()
    print(('finished in %s seconds' % (t1 - t0)))

    print('.. export permutation results')
    perm_file = out_file + '.perm'
    RV = NP.array([LLR0, pv0]).T
    NP.savetxt(perm_file, RV, delimiter='\t', fmt='%.6f %.6e')

    print('.. load test results')
    file_name = os.path.join(resdir, 'test', '*.res')
    files = glob.glob(file_name)
    LLR = []
    for _file in files:
        print(_file)
        _LLR = pd.DataFrame.from_csv(_file, sep='\t')['LLR'].values
        LLR.append(_LLR)
    LLR = NP.concatenate(LLR)

    print('.. calc pvalues')
    pv = c2m.sf(RV_test[:, -1])[:, NP.newaxis]

    print('.. export test results')
    perm_file = out_file + '.test'
    RV_test = NP.hstack([RV_test, pv])
    NP.savetxt(perm_file, RV_test, delimiter='\t',
               fmt='%d %d %d %d %d %d %.6e %.6e')

