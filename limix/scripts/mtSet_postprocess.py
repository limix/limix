#! /usr/bin/env python
# Copyright(c) 2014, The mtSet developers (Francesco Paolo Casale, Barbara Rakitsch, Oliver Stegle)
# All rights reserved.

from limix_legacy.mtSet.core.postprocessCore import postprocess
from optparse import OptionParser

def entry_point():

    parser = OptionParser()
    parser.add_option("--resdir", dest='resdir', type=str, default='./')
    parser.add_option("--outfile", dest='outfile', type=str, default=None)
    parser.add_option("--manhattan_plot", dest='manhattan',action="store_true",default=False)
    parser.add_option("--tol", dest='tol', type=float, default=4e-3)
    (options, args) = parser.parse_args()
    postprocess(options)
