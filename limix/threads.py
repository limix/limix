_max_nthreads = None


def set_max_nthreads(nthreads):
    """ Set the maximum number of threads.

    Parameters
    ----------
    nthreads : int
        Maximum number of threads.
    """
    import dask
    from multiprocessing.pool import ThreadPool

    global _max_nthreads

    nthreads = int(nthreads)
    if nthreads < 1:
        raise ValueError("Cannot set number of threads smaller than one.")
    _max_nthreads = nthreads
    dask.config.set(pool=ThreadPool(_max_nthreads))


def get_max_nthreads():
    """ Get the maximum number of threads

    Returns
    -------
    int
        Maximum number of threads.
    """
    from multiprocessing import cpu_count

    global _max_nthreads

    if _max_nthreads is None:
        return cpu_count()
    return _max_nthreads
