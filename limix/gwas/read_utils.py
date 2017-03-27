def read_geno(geno, i1, i2, type):
    if type=='bed':
        return geno[i1:i2].compute().T
    if type=='hdf5':
        return geno[i1:i2,:].T
