from __future__ import absolute_import as _


def set_coord(x, dim, values):
    r""" Assign a new coordinate or subset an existing one. """
    if dim not in x.coords:
        return x.assign_coords(**{dim: list(values)})
    return x.loc[{dim: x.get_index(dim).isin(values)}]


def take(x, indices, dim):
    r""" Subset a data array on an arbitrary dimension. """
    sl = [slice(None)] * x.ndim
    axis = next(i for i, d in enumerate(x.dims) if d == dim)
    sl[axis] = indices
    return x[tuple(sl)]
