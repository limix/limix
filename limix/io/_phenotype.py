def fetch_bimbam_phenotype(filepath, verbose):
    from .bimbam import read_phenotype
    from xarray import DataArray

    X = read_phenotype(filepath, verbose=verbose)
    return DataArray(X, dims=["sample", "trait"])


def fetch_csv_phenotype(filepath, verbose):
    from .csv import read

    return read(filepath, verbose=verbose)


_dispatch = {"bimbam-pheno": fetch_bimbam_phenotype, "csv": fetch_csv_phenotype}


def fetch_phenotype(fetch_spec, verbose=True):
    from ast import literal_eval

    filetype = fetch_spec["filetype"]
    X = _dispatch[filetype](fetch_spec["filepath"], verbose=verbose)

    spec = fetch_spec["matrix_spec"]
    if spec is not None:
        spec = spec.replace("trait", "trait=")
        spec = dict(
            (k, literal_eval(v))
            for k, v in (pair.split("=") for pair in spec.split())
        )
        X = X.sel(**spec)

    return X
