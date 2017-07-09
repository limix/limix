def entry_point():
    import os
    import time
    from optparse import OptionParser

    import numpy as np
    import pandas as pd
    import scipy as sp

    from ..data import BedReader
    from ..iSet.iset import fit_iSet
    from ..util import unique_variants as f_uni_variants

    parser = OptionParser()
    parser.add_option("--bfile", dest='bfile', type=str, default=None)
    # parser.add_option("--cfile", dest='cfile', type=str, default=None)
    parser.add_option("--pfile", dest='pfile', type=str, default=None)
    parser.add_option("--wfile", dest='wfile', type=str, default=None)
    parser.add_option("--ffile", dest='ffile', type=str, default=None)
    parser.add_option("--ifile", dest='ifile', type=str, default=None)
    parser.add_option("--resdir", dest='resdir', type=str, default='./')

    # start window, end window and permutations
    parser.add_option("--n_perms", type=int, default=10)
    parser.add_option("--start_wnd", dest='i0', type=int, default=None)
    parser.add_option("--end_wnd", dest='i1', type=int, default=None)
    parser.add_option("--factr", dest='factr', type=float, default=1e7)

    parser.add_option(
        "--unique_variants",
        action="store_true",
        dest='unique_variants',
        default=False)

    parser.add_option(
        "--standardize",
        action="store_true",
        dest='standardize',
        default=False)

    (options, args) = parser.parse_args()

    print('importing data')
    F = sp.loadtxt(options.ffile + '.fe')
    Y = sp.loadtxt(options.pfile + '.phe')
    if len(Y.shape) == 1:
        Y = Y[:, sp.newaxis]

    sets = pd.DataFrame.from_csv(
        options.wfile + '.wnd', sep='\t', index_col=None)

    reader = BedReader(options.bfile)

    i0 = 1 if options.i0 is None else options.i0
    i1 = sets.shape[0] if options.i1 is None else options.i1

    df = pd.DataFrame()
    df0 = pd.DataFrame()

    if options.ifile is None:
        Ie = None
    else:
        Ie = sp.loadtxt(options.ifile + '.ind').flatten() == 1

    res_dir = options.resdir

    if not os.path.exists(res_dir):
        os.makedirs(res_dir)

    n_digits = len(str(sets.shape[0]))
    fname = str(i0).zfill(n_digits)
    fname += '_' + str(i1).zfill(n_digits)
    resfile = os.path.join(res_dir, fname)

    for wnd_i in range(i0, i1):
        t0 = time.time()

        _set = sets.ix[wnd_i]
        print('.. set %d: %s' % (wnd_i, _set['setid']))

        Xr = reader.getGenotypes(
            pos_start=_set['start'],
            pos_end=_set['end'],
            chrom=_set['chrom'],
            impute=True)

        if options.unique_variants:
            Xr = f_uni_variants(Xr)

        if options.standardize:
            Xr -= Xr.mean(0)
            Xr /= Xr.std(0)
        else:
            # encoding minor as 0
            p = 0.5 * Xr.mean(0)
            Xr[:, p > 0.5] = 2 - Xr[:, p > 0.5]

        Xr /= np.sqrt(Xr.shape[1])
        _df, _df0 = fit_iSet(Y, F=F, Xr=Xr, Ie=Ie, n_nulls=10)
        df = df.append(_df)
        df0 = df0.append(_df0)
        print('Elapsed:', time.time() - t0)

    df.to_csv(resfile + '.iSet.real')
    df0.to_csv(resfile + '.iSet.perm')
