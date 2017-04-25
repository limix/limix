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

    resdir = options.resdir
    out_file = options.outfile
    tol = options.tol

    print('.. load permutation results')
    file_name = os.path.join(resdir, 'perm*', '*.res')
    files = glob.glob(file_name)
    df0 = []
    for _file in files:
        print(_file)
        _df = pd.DataFrame.from_csv(_file, sep='\t', index_col=None)
        df0.append(_df)
    df0 = pd.concat(df0)

    print('.. fit test statistics')
    t0 = time.time()
    c2m = Chi2mixture(tol=4e-3)
    c2m.estimate_chi2mixture(df0['LLR'].values)
    t1 = time.time()
    print(('finished in %s seconds' % (t1 - t0)))

    print('.. load test results')
    file_name = os.path.join(resdir, 'test', '*.res')
    files = glob.glob(file_name)

    df = []
    for _file in files:
        print(_file)
        _df = pd.DataFrame.from_csv(_file, sep='\t', index_col=None)
        df.append(_df)
    df = pd.concat(df)

    print('.. calc pvalues')
    df0['pv'] = c2m.sf(df0['LLR'].values)
    df['pv'] = c2m.sf(df['LLR'].values)[:, NP.newaxis]

    print('.. export test results')
    perm_file = out_file + '.test'
    df.to_csv(perm_file, sep='\t', index=False)

    print('.. export permutation results')
    perm_file = out_file + '.perm'
    df0.to_csv(perm_file, sep='\t', index=False)
