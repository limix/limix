def fetch_bed_genotype(filepath, verbose=True):
    from .plink import read

    candidates, samples, G = read(filepath, verbose=verbose)

    G = _dataarray_upcast(G)
    G = G.rename(dim_0="candidate", dim_1="sample").T

    for colname in samples.columns:
        G.coords[colname] = ("sample", samples[colname].values)

    G.coords[samples.index.name] = ("sample", samples.index)

    for colname in candidates.columns:
        G.coords[colname] = ("candidate", candidates[colname].values)

    G.coords[candidates.index.name] = ("candidate", candidates.index)

    return G


_dispatch = {"bed": fetch_bed_genotype}


def fetch_genotype(fetch_spec, verbose=True):
    from xarray import DataArray

    filetype = fetch_spec["filetype"]

    spec = fetch_spec["matrix_spec"]
    dims = {d: spec[d] for d in ["row", "col"] if d in spec}

    X = _dispatch[filetype](fetch_spec["filepath"], verbose=verbose)
    X = DataArray(X, dims=_read_dims_into(dims, "sample", "candidate"))

    if len(spec["sel"]) > 0:
        X = X.sel(**spec["sel"])

    X.name = "genotype"

    return X


def _read_dims_into(dims, row, col):
    default = {"row": row, "col": col}
    for dim in default.keys():
        if dim in dims:
            default[dim] = dims[dim]
    return [default["row"], default["col"]]
