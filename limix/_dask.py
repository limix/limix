import dask.array as da
from numpy import asarray, any, isnan


def _get_shape_helper(a):
    s = asarray(a.shape, dtype=int)
    return s[len(s) * (None,) + (slice(None),)]


def _get_all_chunk_shapes(a):
    return a.map_blocks(
        _get_shape_helper,
        dtype=int,
        chunks=tuple(len(c) * (1,) for c in a.chunks) + ((a.ndim,),),
        new_axis=a.ndim,
    )


def _get_chunks(a):
    cs = _get_all_chunk_shapes(a)

    c = []
    for i in range(a.ndim):
        s = a.ndim * [0] + [i]
        s[i] = slice(None)
        s = tuple(s)

        c.append(tuple(cs[s]))

    return tuple(c)


def array_shape_reveal(a):
    if any(isnan(a.shape)):
        # Rebuild Dask Array with known chunks
        return da.Array(a.__dask_graph__(), a.name, _get_chunks(a), a.dtype)
    return a
