from numpy import asarray as npy_asarray


def asarray(X):
    import dask.array as da

    if not isinstance(X, da.Array):
        X = npy_asarray(X, float)
    return X
