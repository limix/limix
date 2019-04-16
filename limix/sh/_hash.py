def filehash(filepath):
    """
    Compute sha256 from a given file.

    Parameters
    ----------
    filepath : str
        File path.

    Returns
    -------
    sha256 : str
        Sha256 of a given file.
    """
    import hashlib

    BUF_SIZE = 65536
    sha256 = hashlib.sha256()

    with open(filepath, "rb") as f:
        while True:
            data = f.read(BUF_SIZE)
            if not data:
                break
            sha256.update(data)

    return sha256.hexdigest()
