from .._data.conf import is_data_name, get_data_dims


def fetch(data_name, fetch_spec, verbose=True):
    from .._data import to_dataarray
    from .._bits.xarray import hint_aware_sel

    if not is_data_name(data_name):
        raise ValueError("`{}` is not a valid data name.".format(data_name))

    filetype = fetch_spec["filetype"]

    spec = fetch_spec["matrix_spec"]
    dims = {d: spec[d] for d in ["row", "col"] if d in spec}

    X = _dispatch[data_name][filetype](fetch_spec["filepath"], verbose=verbose)
    X = to_dataarray(X)
    X = _read_dims_into(X, dims)

    if len(spec["sel"]) > 0:
        X = hint_aware_sel(X, **spec["sel"])

    if X.name is None:
        X.name = data_name

    return X


def _fetch_bed_genotype(filepath, verbose=True):
    from .plink import read
    from xarray import DataArray

    candidates, samples, G = read(filepath, verbose=verbose)

    G = DataArray(G.T, dims=get_data_dims("genotype"))

    for colname in samples.columns:
        G.coords[colname] = ("sample", samples[colname].values)

    G.coords[samples.index.name] = ("sample", samples.index)

    for colname in candidates.columns:
        G.coords[colname] = ("candidate", candidates[colname].values)

    G.coords[candidates.index.name] = ("candidate", candidates.index)

    return G


def _fetch_bimbam_phenotype(filepath, verbose):
    from .bimbam import read_phenotype

    return read_phenotype(filepath, verbose=verbose)


def _fetch_csv_phenotype(filepath, verbose):
    from .csv import read

    return read(filepath, verbose=verbose)


def _read_dims_into(X, dims):
    rc = {"row": 0, "col": 1}
    for dim in dims:
        if rc[dim] != next(i for i in range(len(X.dims)) if X.dims[i] == dim):
            X = X.T
    return X


_dispatch = {
    "genotype": {"bed": _fetch_bed_genotype},
    "trait": {"bimbam-pheno": _fetch_bimbam_phenotype, "csv": _fetch_csv_phenotype},
}
