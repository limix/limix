import os
import pdb
import sys
import csv
import numpy as np
from optparse import OptionParser
import time
import limix
from .read_utils import readNullModelFile
from .read_utils import readWindowsFile
from .read_utils import readCovarianceMatrixFile
from .read_utils import readCovariatesFile
from .read_utils import readPhenoFile
from . import plink_reader
import scipy as sp
import warnings

def scan(bfile,Y,cov,null,wnds,minSnps,i0,i1,perm_i,resfile,F,colCovarType_r='lowrank',rank_r=1,factr=1e7):

    if perm_i is not None:
        print(('Generating permutation (permutation %d)'%perm_i))
        np.random.seed(perm_i)
        perm = np.random.permutation(Y.shape[0])

    mtSet = limix.MTSet(Y=Y, S_R=cov['eval'], U_R=cov['evec'], F=F, rank=rank_r)
    mtSet.setNull(null)
    bim = plink_reader.readBIM(bfile,usecols=(0,1,2,3))
    fam = plink_reader.readFAM(bfile,usecols=(0,1))

    print('fitting model')
    wnd_file = csv.writer(open(resfile,'w'),delimiter='\t')
    for wnd_i in range(i0,i1):
        print(('.. window %d - (%d, %d-%d) - %d snps'%(wnd_i,int(wnds[wnd_i,1]),int(wnds[wnd_i,2]),int(wnds[wnd_i,3]),int(wnds[wnd_i,-1]))))
        if int(wnds[wnd_i,-1])<minSnps:
            print('SKIPPED: number of snps lower than minSnps')
            continue
        #RV = bed.read(PositionRange(int(wnds[wnd_i,-2]),int(wnds[wnd_i,-1])))
        RV = plink_reader.readBED(bfile, useMAFencoding=True, blocksize = 1, start = int(wnds[wnd_i,4]), nSNPs = int(wnds[wnd_i,5]), order  = 'F',standardizeSNPs=False,ipos = 2,bim=bim,fam=fam)

        Xr = RV['snps']
        if perm_i is not None:
            Xr = Xr[perm,:]

        Xr = np.ascontiguousarray(Xr)
        rv = mtSet.optimize(Xr,factr=factr)
        line = np.concatenate([wnds[wnd_i,:],rv['LLR']])
        wnd_file.writerow(line)
    pass

def analyze(options):

    # load data
    print('import data')
    if options.cfile is None:
        cov = {'eval':None,'evec':None}
        warnings.warn('warning: cfile not specifed, a one variance compoenent model will be considered')
    else:
        cov = readCovarianceMatrixFile(options.cfile,readCov=False)
    Y = readPhenoFile(options.pfile,idx=options.trait_idx)
    null = readNullModelFile(options.nfile)
    wnds = readWindowsFile(options.wfile)

    F = None
    if options.ffile:
        F = readCovariatesFile(options.ffile)
        #null['params_mean'] = sp.loadtxt(options.nfile + '.f0')


    if F is not None: assert Y.shape[0]==F.shape[0], 'dimensions mismatch'


    if options.i0 is None: options.i0 = 1
    if options.i1 is None: options.i1 = wnds.shape[0]

    # name of output file
    if options.perm_i is not None:
        res_dir = os.path.join(options.resdir,'perm%d'%options.perm_i)
    else:
        res_dir = os.path.join(options.resdir,'test')
    if not os.path.exists(res_dir):
        os.makedirs(res_dir)
    n_digits = len(str(wnds.shape[0]))
    fname = str(options.i0).zfill(n_digits)
    fname+= '_'+str(options.i1).zfill(n_digits)+'.res'
    resfile = os.path.join(res_dir,fname)

    # analysis
    t0 = time.time()
    scan(options.bfile,Y,cov,null,wnds,options.minSnps,options.i0,options.i1,options.perm_i,resfile,F,options.colCovarType_r,options.rank_r,options.factr)
    t1 = time.time()
    print(('... finished in %s seconds'%(t1-t0)))
