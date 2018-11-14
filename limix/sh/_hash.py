from __future__ import division


# # TODO: document those functions
# def array_hash(X):
#     writeable = X.flags.writeable
#     X.flags.writeable = False
#     h = hash(X.tobytes())
#     X.flags.writeable = writeable
#     return h


def filehash(filepath):
    r"""Compute sha256 from a given file."""
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
