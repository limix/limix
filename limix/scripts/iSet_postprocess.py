#! /usr/bin/env python
# Copyright(c) 2014, The mtSet developers (Francesco Paolo Casale, Barbara Rakitsch, Oliver Stegle)
# All rights reserved.

from optparse import OptionParser
from limix_legacy.mtSet.core.iset_utils import calc_emp_pv_eff
import pandas as pd
import glob
import os
import time
import sys

def entry_point():

    parser = OptionParser()
    parser.add_option("--resdir", dest='resdir', type=str, default='./')
    parser.add_option("--outfile", dest='outfile', type=str, default=None)
    #parser.add_option("--manhattan_plot", dest='manhattan',action="store_true",default=False)
    parser.add_option("--tol", dest='tol', type=float, default=4e-3)
    (options, args) = parser.parse_args()


    resdir = options.resdir
    out_file = options.outfile
    tol = options.tol

    print('.. load permutation results')
    file_name = os.path.join(resdir, '*.iSet.perm')
    files = glob.glob(file_name)
    df0 = pd.DataFrame()
    for _file in files:
        print(_file)
        df0 = df0.append(pd.read_csv(_file, index_col=0))

    print('.. load real results')
    file_name = os.path.join(resdir, '*.iSet.real')
    files = glob.glob(file_name)
    df = pd.DataFrame()
    for _file in files:
        print(_file)
        df = df.append(pd.read_csv(_file, index_col=0))

    #calculate P values for the three tests
    for test in ['mtSet', 'iSet', 'iSet-het']:
        df[test+' pv'] = calc_emp_pv_eff(df[test+' LLR'].values,
                                         df0[test+' LLR0'].values)

    print(('.. saving %s' % out_file+'.res'))
    df.to_csv(out_file+'.res')

    if 0:
        if options.manhattan:
            import limix_legacy.utils.plot as plot

            if not os.path.exists(options.outfile):
                os.makedirs(options.outfile)

            def plot_manhattan(pv, out_file):
                import matplotlib.pylab as PLT
                import scipy as SP
                posCum = SP.arange(pv.shape[0])
                idx=~SP.isnan(pv)
                plot.plot_manhattan(posCum[idx],pv[idx],alphaNS=1.0,alphaS=1.0)
                PLT.savefig(out_file)

            for test in ['mtSet', 'iSet', 'iSet-het']:
                out_file = os.path.join(options.outfile,
                                        'iSet.%s_pv.manhattan.png'\
                                        % (test,))
                print((".. saving " + out_file))
                plot_manhattan(df['%s pv' % test].values, out_file)
