from ._conf import CONF


def get_dims_from_data_name(name):
    return CONF["data_dims"][name]


def is_dim_name(name):
    return name in CONF["dim_names"]


def get_dim_axis(name):
    return CONF["dim_axis"][name]


def is_dim_hint(dim):
    return dim[0] == "_" and is_dim_name(dim[1:])


def dim_name_to_hint(dim):
    if not (is_dim_name(dim) or is_dim_hint(dim)):
        raise ValueError("`{}` is not a dimensional name/hint.")
    if is_dim_name(dim):
        return "_" + dim
    return dim


def dim_hint_to_name(dim):
    if not (is_dim_name(dim) or is_dim_hint(dim)):
        raise ValueError("`{}` is not a dimensional name/hint.")
    if is_dim_hint(dim):
        return dim[1:]
    return dim
