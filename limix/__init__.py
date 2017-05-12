r"""
*************
limix package
*************

A flexible and fast mixed model toolbox.

"""

from __future__ import absolute_import as _absolute_import

import limix_core as core
from pkg_resources import DistributionNotFound as _DistributionNotFound
from pkg_resources import get_distribution as _get_distribution

from . import (heritability, io, iSet, mtSet, plot, qtl, scripts, stats, util,
               varDecomp)
from .mtSet import MTSet

try:
    __version__ = _get_distribution('limix').version
except _DistributionNotFound:
    __version__ = 'unknown'


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
    'test', 'core', 'io', 'plot', 'qtl', 'stats', 'util', 'varDecomp', 'mtSet',
    'iSet', 'scripts', 'MTSet', 'heritability'
]
