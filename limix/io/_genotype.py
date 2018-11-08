def fetch_bed_genotype(filepath):
    from .plink import read
    from limix._dataset import _dataarray_upcast

    candidate, samples, G = read(filepath)

    G = _dataarray_upcast(G)
    G = G.rename(dim_0="candidate", dim_1="sample").T

    for colname in samples.columns:
        G.coords[colname] = ("sample", samples[colname].values)

    for colname in candidate.columns:
        G.coords[colname] = ("candidate", candidate[colname].values)

    return G


_dispatch = {"bed": fetch_bed_genotype}


def fetch_genotype(fetch_spec):
    from ast import literal_eval

    filetype = fetch_spec["filetype"]
    X = _dispatch[filetype](fetch_spec["filepath"])

    spec = fetch_spec["matrix_spec"]
    if spec is not None:
        spec = spec.replace("trait", "trait=").replace("genotype", "genotype=")
        spec = dict(
            (k, literal_eval(v)) for k, v in (pair.split("=") for pair in spec.split())
        )
        X = X.sel(**spec)

    return X
