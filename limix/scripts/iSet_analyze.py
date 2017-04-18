#! /usr/bin/env python
# Copyright(c) 2014, The mtSet developers (Francesco Paolo Casale, Barbara Rakitsch, Oliver Stegle)
# All rights reserved.

import time
import sys
import os
from limix_legacy.iSet.iset import fit_iSet
from optparse import OptionParser
import numpy as np
import pandas as pd
import scipy as sp
import csv

from ..mtSet.core.read_utils import readNullModelFile
from ..mtSet.core.read_utils import readWindowsFile
from ..mtSet.core.read_utils import readCovarianceMatrixFile
from ..mtSet.core.read_utils import readCovariatesFile
from ..mtSet.core.read_utils import readPhenoFile
from ..mtSet.core import plink_reader

def entry_point():
    parser = OptionParser()
    parser.add_option("--bfile", dest='bfile', type=str, default=None)
    #parser.add_option("--cfile", dest='cfile', type=str, default=None)
    parser.add_option("--pfile", dest='pfile', type=str, default=None)
    parser.add_option("--wfile", dest='wfile', type=str, default=None)
    parser.add_option("--ffile", dest='ffile', type=str, default=None)
    parser.add_option("--ifile", dest='ifile', type=str, default=None)
    parser.add_option("--resdir", dest='resdir', type=str, default='./')

    # start window, end window and permutations
    parser.add_option("--minSnps", dest='minSnps', type=int, default=4)
    parser.add_option("--n_perms", type=int, default=10)
    parser.add_option("--start_wnd", dest='i0', type=int, default=None)
    parser.add_option("--end_wnd", dest='i1', type=int, default=None)
    parser.add_option("--factr", dest='factr', type=float, default=1e7)

    (options, args) = parser.parse_args()

    print('importing data')
    F = sp.loadtxt(options.ffile+'.fe')
    Y = sp.loadtxt(options.pfile+'.phe')
    if len(Y.shape)==1: Y = Y[:,sp.newaxis]

    wnds = readWindowsFile(options.wfile)

    bim = plink_reader.readBIM(options.bfile,usecols=(0,1,2,3))
    fam = plink_reader.readFAM(options.bfile,usecols=(0,1))
    chrom = bim[:, 0].astype(float)
    pos = bim[:, -1].astype(float)

    i0 = 1 if options.i0 is None else options.i0
    i1 = wnds.shape[0] if options.i1 is None else options.i1

    df = pd.DataFrame()
    df0 = pd.DataFrame()

    if options.ifile is None:
        Ie = None
    else:
        Ie = sp.loadtxt(options.ifile+'.ind').flatten()==1

    res_dir = options.resdir

    if not os.path.exists(res_dir):
        os.makedirs(res_dir)

    n_digits = len(str(wnds.shape[0]))
    fname = str(i0).zfill(n_digits)
    fname+= '_'+str(i1).zfill(n_digits)
    resfile = os.path.join(res_dir, fname)

    for wnd_i in range(i0,i1):
        t0 = time.time()
        print(('.. window %d - (%d, %d-%d) - %d snps'%(wnd_i,int(wnds[wnd_i,1]),int(wnds[wnd_i,2]),int(wnds[wnd_i,3]),int(wnds[wnd_i,-1]))))
        Xr = plink_reader.readBED(options.bfile, useMAFencoding=True, start = int(wnds[wnd_i,4]), nSNPs = int(wnds[wnd_i,5]), bim=bim , fam=fam)['snps']
        Xr = np.ascontiguousarray(Xr)
        xr = sp.dot(sp.rand(Xr.shape[0]), Xr)
        idxs_u = sp.sort(sp.unique(xr, return_index=True)[1])
        if idxs_u.shape[0]<options.minSnps:
            print('SKIPPED: number of snps lower than minSnps')
            continue
        Xr = Xr[:,idxs_u]
        Xr-= Xr.mean(0)
        Xr/= Xr.std(0)
        Xr/= np.sqrt(Xr.shape[1])
        _df, _df0 = fit_iSet(Y, F=F, Xr=Xr, Ie=Ie, n_nulls=10)
        df  = df.append(_df)
        df0 = df0.append(_df0)
        print('Elapsed:', time.time()-t0)

    df.to_csv(resfile + '.iSet.real')
    df0.to_csv(resfile + '.iSet.perm')
