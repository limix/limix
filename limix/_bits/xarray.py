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


def in_coords_dim(arr, k):
    return k in arr.coords or k in arr.dims


def hint_aware_sel(x, **kwargs):
    from .._data import is_dim_hint, is_dim_name, dim_name_to_hint, dim_hint_to_name

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


def query(data, expr):
    from io import StringIO
    from tokenize import generate_tokens, OP, NAME

    tokens = list(generate_tokens(StringIO(expr).readline))

    final_expr = ""
    last = None
    for t in tokens:
        if t.type == NAME:

            is_boolean = last is not None
            is_boolean = is_boolean and not (last.type == OP and _is_comp(last.string))
            is_boolean = is_boolean and _is_boolean(t.string)
            if is_boolean:
                final_expr += _cast_boolean(t.string)
            else:
                final_expr += 'data["{}"]'.format(t.string)
        elif t.type == OP and _is_comp(t.string):
            final_expr += " {} ".format(t.string)
        else:
            final_expr += t.string
        last = t

    return eval("data.where(" + final_expr + ", drop=True)")


def _is_comp(v):
    return v in set(["<", ">", "<=", ">=", "==", "!="])


def _is_boolean(v):
    return v.lower() in set(["and", "or", "not"])


def _cast_boolean(v):
    d = {"and": " & ", "or": " | ", "not": "~"}
    return d[v.lower()]
