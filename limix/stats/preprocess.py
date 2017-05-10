from __future__ import division

import scipy.spatial
from numpy import (ascontiguousarray, double, einsum, logical_not, newaxis,
                   sqrt, zeros)
from scipy.spatial import _distance_wrap
from tqdm import tqdm
from joblib import Parallel, delayed


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

    import dask
    from chest import Chest

    import tempfile
    import shutil

    #dirpath = tempfile.mkdtemp()
    #cache = Chest(path=dirpath, available_memory=8e8)
    #with dask.set_options(cache=cache):
    #import pdb; pdb.set_trace()

    with tqdm(total=n, desc='Indep. pairwise',
              disable=not verbose) as pbar:

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
                if len(delayeds) == 12:
                    Parallel(n_jobs=12)(delayeds)
                    pbar.update(len(delayeds))
                    delayeds = []


            #Parallel(backend="threading")(delayeds)
            #import pdb; pdb.set_trace()
            Parallel(n_jobs=12)(delayeds)
            pbar.update(len(delayeds))

    #shutil.rmtree(dirpath)
    # import atexit
    # atexit.register(lambda: shutil.rmtree(dirpath))

    return logical_not(excls)
