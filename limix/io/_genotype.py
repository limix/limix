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


def fetch_genotype(filepath, filetype):
    return _dispatch[filetype](filepath)
