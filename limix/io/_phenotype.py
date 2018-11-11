def fetch_bimbam_phenotype(filepath, dims, verbose):
    from .bimbam import read_phenotype
    from xarray import DataArray

    X = read_phenotype(filepath, verbose=verbose)
    X = DataArray(X, dims=_read_dims_into(dims, "sample", "trait"))
    return X


def fetch_csv_phenotype(filepath, dims, verbose):
    from .csv import read
    from xarray import DataArray

    X = read(filepath, verbose=verbose)
    X = DataArray(X, dims=_read_dims_into(dims, "sample", "trait"))
    return X


_dispatch = {"bimbam-pheno": fetch_bimbam_phenotype, "csv": fetch_csv_phenotype}


def fetch_phenotype(fetch_spec, verbose=True):
    filetype = fetch_spec["filetype"]

    spec = fetch_spec["matrix_spec"]
    dims = {d: spec[d] for d in ["row", "col"] if d in spec}

    X = _dispatch[filetype](fetch_spec["filepath"], dims, verbose=verbose)

    if len(spec["sel"]) > 0:
        X = X.sel(**spec["sel"])

    X.name = "phenotype"

    return X


def _read_dims_into(dims, row, col):
    default = {"row": row, "col": col}
    for dim in default.keys():
        if dim in dims:
            default[dim] = dims[dim]
    return [default["row"], default["col"]]
