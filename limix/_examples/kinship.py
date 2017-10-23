import limix
from os.path import join
from numpy import load


def get_kinship(verbose=True):
    url = "http://rest.s3for.me/limix/1000G_kinship.npy"
    with limix.util.TmpDir() as folder:
        limix.util.download(url, folder, verbose=verbose)
        return load(join(folder, "1000G_kinship.npy"))
