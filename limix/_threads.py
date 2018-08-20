from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool

_max_nthreads = None


def set_max_nthreads(nthreads):
    r"""Set the maximum number of threads.
    Parameters
    ----------
    nthreads : int
        Maximum number of threads.
    """
    import dask

    nthreads = int(nthreads)
    if nthreads < 1:
        raise ValueError("Cannot set number of threads smaller than one.")
    _max_nthreads = nthreads
    dask.config.set(pool=ThreadPool(_max_nthreads))


def get_max_nthreads():
    r"""Get the maximum number of threads."""
    if _max_nthreads is None:
        return cpu_count()
    return _max_nthreads
