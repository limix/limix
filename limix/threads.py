from multiprocessing.pool import ThreadPool

import dask
from joblib import cpu_count

_max_nthreads = max(1, cpu_count())
dask.config.set(pool=ThreadPool(_max_nthreads))


def set_max_nthreads(nthreads):
    r"""Set the maximum number of threads.

    Parameters
    ----------
    nthreads : int
        Maximum number of threads.
    """
    global _max_nthreads
    nthreads = int(nthreads)
    if nthreads < 1:
        raise ValueError("Cannot set number of threads smaller than one.")
    _max_nthreads = nthreads

    dask.config.set(pool=ThreadPool(nthreads))


def get_max_nthreads():
    r"""Get the maximum number of threads."""
    return _max_nthreads
