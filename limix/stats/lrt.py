from __future__ import division

from numpy import abs as npy_abs
from numpy import asarray, sqrt
from scipy.stats import chi2


def lrt_pvalues(null_lml, alt_lmls):
    from scipy.stats import chi2
    lrs = -2 * null_lml + 2 * asarray(alt_lmls)
    return chi2(df=1).sf(lrs)


def effsizes_se(effsizes, pvalues):
    return npy_abs(effsizes) / sqrt(chi2(1).isf(pvalues))
