def fetch_dosage(filepath, verbose=True):
    import bgen_reader

    bgen = bgen_reader.read_bgen(filepath, verbose=verbose)
    dosage = bgen_reader.convert_to_dosage(bgen["genotype"], verbose=verbose)
    return dosage.T
