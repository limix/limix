from __future__ import division


def indep_pairwise(X, window_size, step_size, threshold, verbose=True):
    """
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
    ok : ndarray
        Boolean array defining independent variants

    Example
    -------

    .. doctest::

        >>> from numpy.random import RandomState
        >>> from limix.qc import indep_pairwise
        >>>
        >>> random = RandomState(0)
        >>> X = random.randn(10, 20)
        >>>
        >>> indep_pairwise(X, 4, 2, 0.5, verbose=False)
        array([ True,  True, False,  True,  True,  True,  True,  True,  True,
                True,  True,  True,  True,  True,  True,  True,  True,  True,
                True,  True])
    """
    from joblib import Parallel, delayed
    from tqdm import tqdm
    from ..threads import get_max_nthreads
    from numpy import ascontiguousarray, logical_not, zeros

    left = 0
    excls = zeros(X.shape[1], dtype=bool)

    if step_size > window_size:
        raise ValueError("Window size has to be smaller than step size.")

    n = (X.shape[1] + step_size) // step_size

    steps = list(range(n))
    cc = get_max_nthreads()

    with tqdm(total=n, desc="Indep. pairwise", disable=not verbose) as pbar:

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

                delayeds.append(delayed(_func)(x, excls[left:right], threshold))
                if len(delayeds) == cc:
                    Parallel(n_jobs=min(len(delayeds), cc), backend="threading")(
                        delayeds
                    )
                    pbar.update(len(delayeds))
                    delayeds = []

            if len(delayeds) == 0:
                continue

            Parallel(n_jobs=min(len(delayeds), cc), backend="threading")(delayeds)
            pbar.update(len(delayeds))

    return logical_not(excls)


def _row_norms(X):
    from numpy import double, einsum, sqrt

    norms = einsum("ij,ij->i", X, X, dtype=double)
    return sqrt(norms, out=norms)


def _sq_pearson(X):
    from numpy import ascontiguousarray, double, newaxis, zeros
    from scipy.spatial import _distance_wrap

    m = X.shape[0]
    dm = zeros((m * (m - 1)) // 2, dtype=double)

    X2 = X - X.mean(1)[:, newaxis]
    X2 = ascontiguousarray(X2)
    if hasattr(_distance_wrap, "pdist_cosine_wrap"):
        norms = _row_norms(X2)

    X2 = X - X.mean(axis=1, keepdims=True)

    if hasattr(_distance_wrap, "pdist_cosine_wrap"):
        _distance_wrap.pdist_cosine_wrap(X2, dm, norms)
    else:
        _distance_wrap.pdist_cosine_double_wrap(X2, dm)

    return (-dm + 1) ** 2


def _pdist_threshold(mark, dist, thr):
    mark[:] = False
    size = len(mark)

    m = 0
    for i in range(0, size - 1):
        if mark[i]:
            m += size - (i + 1)
            continue

        for j in range(i + 1, size):
            if dist[m] > thr:
                mark[j] = True
            m += 1


def _func(x, excls, threshold):
    from numpy import zeros

    dist = _sq_pearson(x)
    e = zeros(x.shape[0], dtype=bool)
    _pdist_threshold(e, dist, threshold)
    excls |= e
