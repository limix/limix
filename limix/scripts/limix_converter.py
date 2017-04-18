#! /usr/bin/env python
# Copyright(c) 2014, The LIMIX developers (Christoph Lippert, Paolo Francesco Casale, Oliver Stegle)
# All rights reserved.
#
# LIMIX is provided under a 2-clause BSD license.
# See LICENSE for the complete license.

import time
import h5py
import pdb
from limix_legacy.deprecated.io.conversion import *




def entry_point():
    infostring = "limix_conveter.py, Copyright(c) 2014, The LIMIX developers\nlast modified: %s" % time.ctime(os.path.getmtime(__file__))
    print (infostring)

    runner = LIMIX_converter(infostring=infostring)
    (options,args) = runner.parse_args()
    runner.run()
