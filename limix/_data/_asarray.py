def asarray(x, target, dims=None):
    from xarray import DataArray
    from limix._bits.dask import is_dataframe as is_dask_dataframe
    from limix._bits.dask import is_array as is_dask_array
    from limix._bits.dask import is_series as is_dask_series
    from limix._bits.dask import array_shape_reveal
    from limix._bits.xarray import is_dataarray
    from ._conf import CONF
    from numpy import issubdtype, integer

    if target not in CONF["targets"]:
        raise ValueError(f"Unknown target name: {target}.")

    import dask.array as da
    import xarray as xr

    if is_dask_dataframe(x) or is_dask_series(x):
        xidx = x.index.compute()
        x = da.asarray(x)
        x = array_shape_reveal(x)
        x0 = xr.DataArray(x)
        x0.coords[x0.dims[0]] = xidx
        if is_dask_dataframe(x):
            x0.coords[x0.dims[1]] = x.columns
        x = x0
    elif is_dask_array(x):
        x = array_shape_reveal(x)
        x = xr.DataArray(x)

    if not is_dataarray(x):
        x = DataArray(x)

    x.name = target

    while x.ndim < 2:
        rdims = set(CONF["data_dims"]["trait"]).intersection(set(x.coords.keys()))
        rdims = rdims - set(x.dims)
        if len(rdims) == 1:
            dim = rdims.pop()
        else:
            dim = "dim_{}".format(x.ndim)
        x = x.expand_dims(dim, x.ndim)

    if isinstance(dims, (tuple, list)):
        dims = {a: n for a, n in enumerate(dims)}
    dims = _numbered_axes(dims)
    if len(set(dims.values())) < len(dims.values()):
        raise ValueError("`dims` must not contain duplicated values.")

    x = x.rename({x.dims[axis]: name for axis, name in dims.items()})
    x = _set_missing_dim(x, CONF["data_dims"][target])
    x = x.transpose(*CONF["data_dims"][target])

    if issubdtype(x.dtype, integer):
        x = x.astype(float)

    for dim in x.dims:
        if x.coords[dim].dtype.kind in {"U", "S"}:
            x.coords[dim].values = x.coords[dim].values.astype(object)

    return x


def _numbered_axes(dims):

    if dims is None:
        return {}

    naxes = {}
    for a, v in dims.items():
        if a == "row":
            naxes[0] = v
        elif a == "col":
            naxes[1] = v
        else:
            naxes[a] = v

    return naxes


def _set_missing_dim(arr, dims):
    unk_dims = set(arr.dims) - set(dims)
    if len(unk_dims) > 1:
        raise ValueError("Too many unknown dimension names.")
    elif len(unk_dims) == 1:
        known_dims = set(dims) - set(arr.dims)
        if len(known_dims) != 1:
            raise ValueError("Can't figure out what is the missing dimension name.")
        arr = arr.rename({unk_dims.pop(): known_dims.pop()})
    return arr
