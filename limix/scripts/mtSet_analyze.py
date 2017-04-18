#! /usr/bin/env python
# Copyright(c) 2014, The mtSet developers (Francesco Paolo Casale, Barbara Rakitsch, Oliver Stegle)
# All rights reserved.

import sys
import os
from limix_legacy.mtSet.core.analyzeCore import analyze
from optparse import OptionParser

def entry_point():
    parser = OptionParser()
    parser.add_option("--bfile", dest='bfile', type=str, default=None)
    parser.add_option("--cfile", dest='cfile', type=str, default=None)
    parser.add_option("--pfile", dest='pfile', type=str, default=None)
    parser.add_option("--nfile", dest='nfile', type=str, default=None)
    parser.add_option("--wfile", dest='wfile', type=str, default=None)
    parser.add_option("--ffile", dest='ffile', type=str, default=None)
    parser.add_option("--resdir", dest='resdir', type=str, default='./')
    parser.add_option("--trait_idx",dest='trait_idx',type=str, default=None)

    parser.add_option("--rank_r",dest='rank_r',type=int, default=1)
    parser.add_option("--colCovarType_r",dest='colCovarType_r',type=str, default='lowrank')

    # start window, end window and permutations
    parser.add_option("--minSnps", dest='minSnps', type=int, default=None)
    parser.add_option("--start_wnd", dest='i0', type=int, default=None)
    parser.add_option("--end_wnd", dest='i1', type=int, default=None)
    parser.add_option("--perm", dest='perm_i', type=int, default=None)
    parser.add_option("--factr", dest='factr', type=float, default=1e7)

    (options, args) = parser.parse_args()

    analyze(options)
