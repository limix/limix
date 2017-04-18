#! /usr/bin/env python
# Copyright(c) 2014, The mtSet developers (Francesco Paolo Casale, Barbara Rakitsch, Oliver Stegle)
# All rights reserved.

from limix_legacy.mtSet.core.simPhenoCore import simPheno
from optparse import OptionParser
import scipy as SP

def entry_point():

    parser = OptionParser()
    parser.add_option("--bfile",     dest='bfile',      type=str, default=None)
    parser.add_option("--cfile",     dest='cfile',      type=str, default=None)
    parser.add_option("--pfile",     dest='pfile',      type=str, default=None)

    parser.add_option("--seed",      dest='seed',       type=int, default=0)
    parser.add_option("--nTraits",   dest='nTraits',    type=int, default=4)
    parser.add_option("--windowSize",dest='windowSize', type=int, default=3e4)

    parser.add_option("--chrom", dest='chrom',  type=int, default=None)
    parser.add_option("--minPos",dest='pos_min',type=int, default=None)
    parser.add_option("--maxPos",dest='pos_max',type=int, default=None)

    parser.add_option("--vTotR",   dest='vTotR', type=float, default=0.05)
    parser.add_option("--nCausalR",dest='nCausalR', type=int,default=8)
    parser.add_option("--pCommonR",dest='pCommonR', type=float,default=0.5)
    parser.add_option("--vTotBg",  dest='vTotBg',type=float,default=0.4)
    parser.add_option("--pHidden", dest='pHidden',type=float,default=0.6)
    parser.add_option("--pCommon", dest='pCommon', type=float,default=0.5)
    (options, args) = parser.parse_args()

    simPheno(options)
