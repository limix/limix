def fetch_bimbam_phenotype(filepath):
    from .bimbam import read_phenotype

    return read_phenotype(filepath)


def fetch_csv_phenotype(filepath):
    from .csv import read

    return read(filepath)


_dispatch = {"bimbam-pheno": fetch_bimbam_phenotype, "csv": fetch_csv_phenotype}


def fetch_phenotype(fetch_spec):
    filetype = fetch_spec["filetype"]
    df = _dispatch[filetype](fetch_spec["filepath"])

    cols = fetch_spec["matrix_spec"]["cols"]
    if cols != None:
        df = eval("df[" + cols + "]")

    rows = fetch_spec["matrix_spec"]["rows"]
    if rows != None:
        df = eval("df.loc[" + rows + "]")

    return df
