from .._display import session_line


def extract(filepath, verbose=True):
    """
    Extract a compressed file.

    Parameters
    ----------
    filepath : str
        File path.
    verbose : bool, optional
        ``True`` for displaying progress. Defaults to ``True``.
    """
    formats_order = _best_order(filepath)

    err_msgs = []
    with session_line("Extracting {}... ".format(filepath), disable=not verbose):

        for f in formats_order:
            try:
                return f(filepath)
            except Exception as e:
                err_msgs.append(str(e))

        raise RuntimeError(f"Could not extract `{filepath}`.\n" + "\n".join(err_msgs))


def _extract_tar(filepath):
    import tarfile

    tar = tarfile.open(filepath)
    tar.extractall()
    filepath = tar.getnames()[0]
    tar.close()
    return filepath


def _extract_bz2(filepath):
    import bz2
    import os

    filename = os.path.splitext(filepath)[0]

    with open(filepath, "rb") as f:
        o = bz2.decompress(f.read())

    with open(filename, "wb") as f:
        f.write(o)
    return filename


def _extract_zip(filepath, dest="."):
    import zipfile

    with zipfile.ZipFile(filepath, "r") as zip_ref:
        zip_ref.extractall(dest)

    return dest


def _best_order(filepath):
    formats = [_extract_tar, _extract_bz2, _extract_zip]
    formats_order = []
    if filepath.endswith(".zip"):
        formats_order = [_extract_zip]
    elif filepath.endswith(".tar.gz"):
        formats_order = [_extract_tar]
    elif filepath.endswith(".bz2"):
        formats_order = [_extract_bz2]

    for f in formats:
        if f not in formats_order:
            formats_order.append(f)

    return formats_order
