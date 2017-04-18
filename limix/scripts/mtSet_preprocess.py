#! /usr/bin/env python
# Copyright(c) 2014, The mtSet developers (Francesco Paolo Casale, Barbara Rakitsch, Oliver Stegle)
# All rights reserved.

from limix_legacy.mtSet.core.preprocessCore import preprocess
from optparse import OptionParser

def entry_point():
    parser = OptionParser()
    parser.add_option("--bfile", dest='bfile', type=str, default=None)
    parser.add_option("--cfile", dest='cfile', type=str, default=None)
    parser.add_option("--pfile", dest='pfile', type=str, default=None)
    parser.add_option("--nfile", dest='nfile', type=str, default=None)
    parser.add_option("--wfile", dest='wfile', type=str, default=None)
    parser.add_option("--ffile", dest='ffile', type=str, default=None)
    parser.add_option("--trait_idx",dest='trait_idx',type=str, default=None)

    parser.add_option("--compute_covariance", action="store_true",dest="compute_cov",default=False)
    parser.add_option("--compute_PCs", dest="compute_PCs",default=0,type=int)


    parser.add_option("--plink_path", dest='plink_path', type=str, default='plink')
    parser.add_option("--sim_type",   dest='sim_type', type=str, default='RRM')

    parser.add_option("--fit_null", action="store_true",dest="fit_null",default=False)

    parser.add_option("--precompute_windows", action="store_true", dest='precompute_windows', default=False)
    parser.add_option("--window_size", dest='window_size', type=int, default=3e4)
    parser.add_option("--plot_windows", action='store_true', dest='plot_windows',default=False)

    (options, args) = parser.parse_args()

    preprocess(options)
