from __future__ import division

import scipy.spatial
from joblib import Parallel, cpu_count, delayed
from numpy import (asarray, ascontiguousarray, double, einsum, isfinite,
                   logical_not, minimum, newaxis, sqrt, unique, zeros)
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


def func(x, excls, threshold):
    dist = _sq_pearson(x)
    e = zeros(x.shape[0], dtype=bool)
    _pdist_threshold(e, dist, threshold)
    excls |= e


def indep_pairwise(X, window_size, step_size, threshold, verbose=True):
    r"""
    Determine pair-wise independent variants.

    Independent variants are defined via squared Pearson correlations between
    pairs of variants inside a sliding window.

    Parameters
    ----------
    X : array_like
        Sample by variants matrix.
    window_size : int
        Number of variants inside each window.
    step_size : int
        Number of variants the sliding window skips.
    threshold : float
        Squared Pearson correlation threshold for independence.
    verbose : bool
        `True` for progress information; `False` otherwise.

    Returns
    -------
    ok : boolean array defining independent variants

    Examples
    --------
    .. doctest::

        >>> from numpy.random import RandomState
        >>> from limix.stats import indep_pairwise
        >>>
        >>> random = RandomState(0)
        >>> X = random.randn(10, 20)
        >>>
        >>> indep_pairwise(X, 4, 2, 0.5, verbose=False)
        array([ True,  True, False,  True,  True,  True,  True,  True,  True,
                True,  True,  True,  True,  True,  True,  True,  True,  True,
                True,  True], dtype=bool)
    """
    left = 0
    excls = zeros(X.shape[1], dtype=bool)
    excl = zeros(window_size, dtype=bool)

    assert step_size <= window_size

    n = (X.shape[1] + step_size) // step_size

    steps = list(range(n))
    cc = max(1, cpu_count())

    with tqdm(total=n, desc='Indep. pairwise', disable=not verbose) as pbar:

        while len(steps) > 0:
            i = 0
            right = 0
            delayeds = []
            while i < len(steps):

                step = steps[i]
                left = step * step_size
                if left < right:
                    i += 1
                    continue

                del steps[i]
                right = min(left + window_size, X.shape[1])
                x = ascontiguousarray(X[:, left:right].T)

                delayeds.append(delayed(func)(x, excls[left:right], threshold))
                if len(delayeds) == cc:
                    Parallel(
                        n_jobs=min(len(delayeds), cc),
                        backend='threading')(delayeds)
                    pbar.update(len(delayeds))
                    delayeds = []

            if len(delayeds) == 0:
                continue

            Parallel(
                n_jobs=min(len(delayeds), cc), backend='threading')(delayeds)
            pbar.update(len(delayeds))

    return logical_not(excls)


def _check_encoding(X):
    u = unique(X)
    u = u[isfinite(u)]
    if len(u) > 3:
        return False
    return all([i in set([0, 1, 2]) for i in u])


def maf(X):
    r"""Compute minor allele frequencies.

    It assumes that `X` encodes 0, 1, and 2 representing the number
    of alleles.

    Args:
        X (array_like): Genotype matrix.

    Returns:
        array_like: minor allele frequencies.

    Examples
    --------

        .. doctest::

            >>> from numpy.random import RandomState
            >>> from limix.stats import maf
            >>>
            >>> random = RandomState(0)
            >>> X = random.randint(0, 3, size=(100, 10))
            >>>
            >>> print(maf(X))
            [ 0.49   0.49   0.445  0.495  0.5    0.45   0.48   0.48   0.47   0.435]
    """
    ok = _check_encoding(X)
    if not ok:
        raise ValueError("It assumes that X encodes 0, 1, and 2 only.")
    s0 = X.sum(0)
    s0 = s0 / (2 * X.shape[0])
    s1 = 1 - s0
    return minimum(s0, s1)
