def fetch_bimbam_phenotype(filepath, dims, verbose):
    from .bimbam import read_phenotype

    return read_phenotype(filepath, verbose=verbose)


def fetch_csv_phenotype(filepath, dims, verbose):
    from .csv import read

    return read(filepath, verbose=verbose)


_dispatch = {"bimbam-pheno": fetch_bimbam_phenotype, "csv": fetch_csv_phenotype}


def fetch_phenotype(fetch_spec, verbose=True):
    from xarray import DataArray

    filetype = fetch_spec["filetype"]

    spec = fetch_spec["matrix_spec"]
    dims = {d: spec[d] for d in ["row", "col"] if d in spec}

    X = _dispatch[filetype](fetch_spec["filepath"], dims, verbose=verbose)
    X = DataArray(X, dims=_read_dims_into(dims, "sample", "trait"))

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
