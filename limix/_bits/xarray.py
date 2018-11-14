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


def in_coords_dim(arr, k):
    return k in arr.coords or k in arr.dims

def hint_aware_sel(x, **kwargs):
    from .._data.conf import (
        is_dim_hint,
        is_dim_name,
        dim_name_to_hint,
        dim_hint_to_name,
    )

    for k in kwargs.keys():
        if in_coords_dim(x, k):
            continue
        if is_dim_name(k) or is_dim_hint(k):
            if in_coords_dim(x, dim_name_to_hint(k)):
                new_k = dim_name_to_hint(k)
                if new_k not in kwargs:
                    kwargs[new_k] = kwargs[k]
                    del kwargs[k]
            elif in_coords_dim(x, dim_hint_to_name(k)):
                new_k = dim_hint_to_name(k)
                if new_k not in kwargs:
                    kwargs[new_k] = kwargs[k]
                    del kwargs[k]

    return x.sel(**kwargs)
