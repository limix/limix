from __future__ import division

from numpy import isnan


def get_imputer(m):
    def impute(X):
        A = X.copy()

        isn = isnan(A)
        A[:] = 0
        A[isn] = 1

        X[isn] = 0
        X += A * m

        return X

    return impute


def mean_impute(G):
    import dask.array as da

    if isinstance(G, da.Array):
        m = da.nanmean(G, axis=0).compute()
        start = 0

        arrs = []
        for i in range(len(G.chunks[1])):
            end = start + G.chunks[1][i]
            impute = get_imputer(m[start:end])
            arrs.append(G[:, start:end].map_blocks(impute, dtype=float))
            start = end
        G = da.concatenate(arrs, axis=1)
    else:
        raise NotImplementedError

    return G
