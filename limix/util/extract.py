import tarfile


def extract(filepath, verbose=True):
    if verbose:
        print("Extracting {}...".format(filepath))
    tar = tarfile.open(filepath)
    tar.extractall()
    tar.close()
