__all__ = ["set_coord", "take", "query"]


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


def query(x, expr):
    from io import StringIO
    from tokenize import generate_tokens, OP, NAME

    expr = expr.strip()
    if len(expr) == 0:
        return x
    tokens = list(generate_tokens(StringIO(expr).readline))

    final_expr = ""
    last = None
    for t in tokens:
        if t.type == NAME and t.string not in ["isin", "None"]:

            is_boolean = last is not None
            is_boolean = is_boolean and not (last.type == OP and _is_comp(last.string))
            is_boolean = is_boolean and _is_boolean(t.string)
            if is_boolean:
                final_expr += _cast_boolean(t.string)
            else:
                final_expr += 'x["{}"]'.format(t.string)
        elif t.type == OP and _is_comp(t.string):
            final_expr += " {} ".format(t.string)
        else:
            final_expr += t.string
        last = t

    return eval("x.where(" + final_expr + ", drop=True)")


def is_dataarray(a):
    pkg = a.__class__.__module__.split(".")[0]
    name = a.__class__.__name__

    return pkg == "xarray" and name == "DataArray"


def _is_comp(v):
    return v in set(["<", ">", "<=", ">=", "==", "!="])


def _is_boolean(v):
    return v.lower() in set(["and", "or", "not"])


def _cast_boolean(v):
    d = {"and": " & ", "or": " | ", "not": "~"}
    return d[v.lower()]
