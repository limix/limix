def fetch_bimbam_phenotype(filepath):
    from .bimbam import read_phenotype

    return read_phenotype(filepath)


_dispatch = {"bimbam-pheno": fetch_bimbam_phenotype}


def fetch_phenotype(filepath, filetype):
    return _dispatch[filetype](filepath)
