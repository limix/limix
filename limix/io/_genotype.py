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
    filetype = fetch_spec["filetype"]
    df = _dispatch[filetype](fetch_spec["filepath"])

    cols = fetch_spec["matrix_spec"]["cols"]
    if cols != None:
        df = eval("df[" + cols + "]")

    rows = fetch_spec["matrix_spec"]["rows"]
    if rows != None:
        df = eval("df.loc[" + rows + "]")

    return df
