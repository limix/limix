import sys
import os
from os.path import dirname, realpath, join
import tempfile
import shutil

_filenames = [
    "data.csv",
    "data.h5.bz2",
    "example.gen",
    "example.sample",
    "pheno.csv",
    "xarr.hdf5.bz2",
    "ex0/phenotype.gemma",
]


class file_example(object):
    r"""Example."""

    def __init__(self, filenames):
        self._unlist = False
        if not isinstance(filenames, (tuple, list)):
            filenames = [filenames]
            self._unlist = True

        for fn in filenames:
            if fn not in _filenames:
                raise ValueError(
                    "Unrecognized file name {}. Choose one of these: {}".format(
                        fn, _filenames
                    )
                )

        self._dirpath = tempfile.mkdtemp()
        self._filenames = filenames

    def __enter__(self):
        import pkg_resources

        filepaths = [join(self._dirpath, fn) for fn in self._filenames]

        for fn, fp in zip(self._filenames, filepaths):
            if __name__ == "__main__":
                shutil.copy(join(dirname(realpath(__file__)), fn), fp)
            else:
                resource_path = "example/{}".format(fn)
                content = pkg_resources.resource_string(
                    __name__.split(".")[0], resource_path
                )

                _makedirs(dirname(fp))
                with open(fp, "wb") as f:
                    f.write(content)

        if self._unlist:
            return filepaths[0]
        return filepaths

    def __exit__(self, *_):
        try:
            shutil.rmtree(self._dirpath)
        except PermissionError:
            pass


def _makedirs(dirpath):
    if sys.version_info >= (3,):
        os.makedirs(dirpath, exist_ok=True)
    else:
        try:
            os.makedirs(dirpath)
        except OSError:
            pass
