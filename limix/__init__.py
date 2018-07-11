r"""
*************
limix package
*************

A flexible and fast mixed model toolbox.

"""

from __future__ import absolute_import as _absolute_import

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import limix_core as core

from . import (heritability, io, iset, mtset, plot, qtl, scripts, stats, util,
               vardec)

__version__ = '1.0.18'


def test():
    import os
    p = __import__('limix').__path__[0]
    src_path = os.path.abspath(p)
    old_path = os.getcwd()
    os.chdir(src_path)

    try:
        return_code = __import__('pytest').main(['-q', '--doctest-modules'])
    finally:
        os.chdir(old_path)

    if return_code == 0:
        print("Congratulations. All tests have passed!")

    return return_code


__all__ = [
    'test', 'core', 'io', 'plot', 'qtl', 'stats', 'util', 'vardec', 'mtset',
    'iset', 'scripts', 'heritability', '__version__'
]
