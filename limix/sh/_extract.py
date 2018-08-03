import bz2
import os
import tarfile


# TODO: document
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

    filename = os.path.splitext(filepath)[0]

    if os.path.exists(filename):
        if verbose:
            print("File {} already exists.".format(filename))
        return

    with open(filepath, "rb") as f:
        o = bz2.decompress(f.read())

    with open(filename, "wb") as f:
        f.write(o)
