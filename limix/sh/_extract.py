from .._display import session_line


# TODO: document
def extract(filepath, verbose=True):
    import bz2
    import tarfile
    import os

    with session_line("Extracting {}... ".format(filepath), disable=not verbose):

        try:
            tar = tarfile.open(filepath)
            tar.extractall()
            filepath = tar.getnames()[0]
            tar.close()
            return filepath
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

        return filename
