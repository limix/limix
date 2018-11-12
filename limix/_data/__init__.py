
def to_dataarray(x):
    import dask.dataframe as dd
    import dask.array as da
    import xarray as xr

    if x is None:
        return None

    if isinstance(x, (dd.Series, dd.DataFrame)):
        xidx = x.index.compute()
        x = da.asarray(x)
        x = array_shape_reveal(x)
        x0 = xr.DataArray(x)
        x0.coords[x0.dims[0]] = xidx
        if isinstance(x, dd.DataFrame):
            x0.coords[x0.dims[1]] = x.columns
        x = x0

    if not isinstance(x, xr.DataArray):
        x = xr.DataArray(x, encoding={"dtype": "float64"})

    if x.dtype != dtype("float64"):
        x = x.astype("float64")

    if x.ndim < 2:
        x = x.expand_dims("dim_1", 1)

    for dim in x.dims:
        if x.coords[dim].dtype.kind in {"U", "S"}:
            x.coords[dim] = x.coords[dim].values.astype(object)
    return x
