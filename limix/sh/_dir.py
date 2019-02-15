def makedirs(dirpath):
    import os
    import sys

    if sys.version_info >= (3,):
        os.makedirs(dirpath, exist_ok=True)
    else:
        try:
            os.makedirs(dirpath)
        except OSError:
            pass
