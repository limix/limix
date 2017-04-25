from limix.util import sets_from_bim
from limix.util import annotate_sets
from limix.io import read_plink
import time


def definesets(opt):

    assert opt.bfile is not None, 'Please specify a bfile.'

    """ precomputing the windows """
    if opt.wfile is None:
        opt.wfile = os.path.split(
            opt.bfile)[-1] + '.%d' % opt.window_size
        warnings.warn('wfile not specifed, set to %s' % opt.wfile)

    t0 = time.time()

    (bim, fam, G) = read_plink(opt.bfile, verbose=False)

    if opt.sliding_wind:
        print('Precomputing windows')
        sets = sets_from_bim(bim,
                             size=opt.window_size,
                             step=opt.step,
                             chrom=opt.chrom,
                             minSnps=opt.minSnps,
                             maxSnps=opt.maxSnps)

    elif opt.filter_sets:
        assert opt.iwfile is not None, 'Please specify a iwfile.'
        print('Filtering windows')
        pdb.set_trace()
        sets = annotate_sets(sets0,
                             bim,
                             minSnps=opt.minSnps,
                             maxSnps=opt.maxSnps)

    sets.to_csv(opt.wfile + '.wnd', sep='\t', index=False)

    print(('Number of variants:', G.shape[0]))
    print(('Number of windows:', sets.shape[0]))
    t1 = time.time()
    print(('.. finished in %s seconds' % (t1 - t0)))
