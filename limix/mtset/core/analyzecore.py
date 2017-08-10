import os
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
from limix.data import BedReader
from limix.util import unique_variants as f_uni_variants
import scipy as sp
import warnings
import pandas as pd


def scan(
        bfile,
        Y,
        cov,
        null,
        sets,
        i0,
        i1,
        perm_i,
        resfile,
        F,
        colCovarType_r='lowrank',
        rank_r=1,
        factr=1e7,
        unique_variants=False,
        standardize=False):

    if perm_i is not None:
        print(('Generating permutation (permutation %d)' % perm_i))
        np.random.seed(perm_i)
        perm = np.random.permutation(Y.shape[0])

    mtSet = limix.MTSet(
        Y=Y,
        S_R=cov['eval'],
        U_R=cov['evec'],
        F=F,
        rank=rank_r)
    mtSet.setNull(null)

    reader = BedReader(bfile)

    wnd_ids = sp.arange(i0, i1)
    LLR = sp.zeros(sets.shape[0])
    for wnd_i in wnd_ids:

        _set = sets.ix[wnd_i]
        print('.. set %d: %s' % (wnd_i, _set['setid']))

        Xr = reader.getGenotypes(pos_start=_set['start'],
                                 pos_end=_set['end'],
                                 chrom=_set['chrom'],
                                 impute=True)

        if unique_variants:
            Xr = f_uni_variants(Xr)

        if standardize:
            Xr -= Xr.mean(0)
            Xr /= Xr.std(0)
        else:
            # encoding minor as 0
            p = 0.5 * Xr.mean(0)
            Xr[:, p > 0.5] = 2 - Xr[:, p > 0.5]

        if perm_i is not None:
            Xr = Xr[perm, :]

        # multi trait set test fit
        RV = mtSet.optimize(Xr, factr=factr)
        LLR[wnd_i] = RV['LLR'][0]

    # export results
    sets['LLR'] = LLR
    sets.to_csv(resfile, sep='\t', index=False)


def analyze(options):

    # load data
    print('import data')
    if options.cfile is None:
        cov = {'eval': None, 'evec': None}
        warnings.warn(
            'warning: cfile not specifed, a one variance compoenent model will be considered')
    else:
        cov = readCovarianceMatrixFile(options.cfile, readCov=False)
    Y = readPhenoFile(options.pfile, idx=options.trait_idx)
    null = readNullModelFile(options.nfile)

    sets = pd.DataFrame.from_csv(options.wfile + '.wnd',
                                 sep='\t',
                                 index_col=None)

    F = None
    if options.ffile:
        F = readCovariatesFile(options.ffile)
        # null['params_mean'] = sp.loadtxt(options.nfile + '.f0')

    if F is not None:
        assert Y.shape[0] == F.shape[0], 'dimensions mismatch'

    if options.i0 is None:
        options.i0 = 1
    if options.i1 is None:
        options.i1 = sets.shape[0]

    # name of output file
    if options.perm_i is not None:
        res_dir = os.path.join(options.resdir, 'perm%d' % options.perm_i)
    else:
        res_dir = os.path.join(options.resdir, 'test')
    if not os.path.exists(res_dir):
        os.makedirs(res_dir)
    n_digits = len(str(sets.shape[0]))
    fname = str(options.i0).zfill(n_digits)
    fname += '_' + str(options.i1).zfill(n_digits) + '.res'
    resfile = os.path.join(res_dir, fname)

    # analysis
    t0 = time.time()
    scan(
        options.bfile,
        Y,
        cov,
        null,
        sets,
        options.i0,
        options.i1,
        options.perm_i,
        resfile,
        F,
        options.colCovarType_r,
        options.rank_r,
        options.factr,
        options.unique_variants,
        options.standardize)
    t1 = time.time()
    print(('... finished in %s seconds' % (t1 - t0)))
