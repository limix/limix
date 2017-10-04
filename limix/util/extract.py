import bz2
import os
import tarfile


def extract(filepath, verbose=True):
    if verbose:
        print("Extracting {}...".format(filepath))

    try:
        tar = tarfile.open(filepath)
        tar.extractall()
        tar.close()
        return
    except tarfile.ReadError:
        pass

    with open(filepath, 'rb') as f:
        o = bz2.decompress(f.read())

    filename, file_extension = os.path.splitext(filepath)

    with open(filename, 'wb') as f:
        f.write(o)
