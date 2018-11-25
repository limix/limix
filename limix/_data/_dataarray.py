from ._dim import is_dim_name, get_dim_axis


def rename_dims(x, dims):
    r""" Rename dimensions while considering hinted name and default matrix layout. """
    from numpy import atleast_1d

    for i, dim in enumerate(dims):
        if x.dims[i] != dim:
            coords = None
            if dim in x.coords:
                coords = x.coords[dim]
                x = x.drop([dim])
            if x.dims[i].startswith("_") and x.dims[i][1:] == dim:
                if dim in x.coords:
                    x = x.drop([x.dims[i]])
                x = x.rename({x.dims[i]: dim})
                if coords is not None:
                    x = x.assign_coords(**{x.dims[i]: atleast_1d(coords)})
            else:
                x = x.rename({x.dims[i]: dim})
                if coords is not None:
                    x = x.assign_coords(**{dim: atleast_1d(coords)})

    return x


def fix_dim_hint(x):
    r""" Hints start with underscore.

    A DataArray with row-dimension name `_trait` and col-dimension name `_sample`
    hints that traits are row-wise stored and samples are col-wise stored.

    Parameters
    ----------
    x : :class:`xarray.DataArray`
        Bi-dimensional data array.
    """
    hint = {}
    for i, d in enumerate(x.dims):
        if _is_dim_hint(d):
            hint[_hint_to_name(d)] = i

    if len(hint) == 0:
        return x

    for d, v in hint.items():
        if v != get_dim_axis(d):
            x = x.T
            break

    return x


def _is_dim_hint(name):
    return name.startswith("_") and is_dim_name(_hint_to_name(name))


def _hint_to_name(name):
    return name[1:]
