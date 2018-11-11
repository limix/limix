def fetch_bed_genotype(filepath, verbose=True):
    from .plink import read
    from limix._dataset import _dataarray_upcast

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
    filetype = fetch_spec["filetype"]
    spec = fetch_spec["matrix_spec"]

    X = _dispatch[filetype](fetch_spec["filepath"], verbose=verbose)

    if len(spec["sel"]) > 0:
        X = X.sel(**spec["sel"])

    X.name = "genotype"
    return X
