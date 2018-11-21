from __future__ import absolute_import as _


def array_shape_reveal(a):
    r""" Compute shape from lazy arrays.

    Dask arrays might have unknown dimensions due to performance
    reasons. This function compute those dimensions.
    """
    from numpy import any, isnan
    import dask.array as da

    if any(isnan(a.shape)):
        # Rebuild Dask Array with known chunks
        return da.Array(a.__dask_graph__(), a.name, _get_chunks(a), a.dtype)
    return a


def _get_shape_helper(a):
    from numpy import asarray

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
    return tuple(c)
