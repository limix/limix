from __future__ import division

import scipy.spatial
from numpy import (ascontiguousarray, double, einsum, logical_not, newaxis,
                   sqrt, zeros)
from scipy.spatial import _distance_wrap
from tqdm import tqdm


def _row_norms(X):
    norms = einsum('ij,ij->i', X, X, dtype=double)
    return sqrt(norms, out=norms)


def _sq_pearson(X):
    m = X.shape[0]
    dm = zeros((m * (m - 1)) // 2, dtype=double)

    X2 = X - X.mean(1)[:, newaxis]
    X2 = ascontiguousarray(X2)
    norms = _row_norms(X2)
    _distance_wrap.pdist_cosine_wrap(X2, dm, norms)
    return (-dm + 1)**2


def _pdist_threshold(mark, dist, thr):
    mark[:] = False
    size = len(mark)

    l = 0
    for i in range(0, size - 1):
        if mark[i]:
            l += size - (i + 1)
            continue

        for j in range(i + 1, size):
            if dist[l] > thr:
                mark[j] = True
            l += 1


# sample by snp
def indep_pairwise(X, window_size, window_step, threshold, verbose=True):
    left = 0
    excls = zeros(X.shape[1], dtype=bool)
    excl = zeros(window_size, dtype=bool)

    assert window_step <= window_size

    n = (X.shape[1] + window_step) // window_step
    for i in tqdm(range(n), desc='Indep. pairwise', disable=not verbose):

        right = min(left + window_size, X.shape[1])
        x = ascontiguousarray(X[:, left:right].T)

        dist = _sq_pearson(x)

        e = excl[:x.shape[0]]
        _pdist_threshold(e, dist, threshold)

        excls[left:right] |= e

        left += window_step

    return logical_not(excls)
